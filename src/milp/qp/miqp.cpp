#include "milp/qp/miqp.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <queue>
#include <string>

#include "milp/heuristics/fj.hpp"
#include "milp/qp/qp.hpp"
#include "milp/qp/qp_active.hpp"
#include "milp/search/propagation.hpp"

namespace milp {

namespace {

constexpr double kIntTol = 1e-6;   // integrality
constexpr double kFeasTol = 1e-6;  // rows and bounds of an incumbent

struct BoundChange {
  Index col;
  double lb, ub;
};

// Warm start of a node's active-set solve: the parent's final working set and point.
struct Warm {
  WorkingSet ws;
  std::vector<double> x;
  uint64_t state = 0;  // the active-set factorization of ws (if still cached)
  long long bytes() const {
    return static_cast<long long>(sizeof(Warm) + ws.size() * sizeof(int32_t) + x.size() * sizeof(double));
  }
};
using WarmPtr = std::shared_ptr<const Warm>;

// Result of a node relaxation.
struct NodeQp {
  QpStatus status = QpStatus::NumericalError;
  std::vector<double> x;
  double objective = kInf;
  double bound = -kInf;   // lower bound on the node's QP (the safe bound when finite)
  std::vector<double> d;  // reduced costs of the safe bound (empty if none)
  WarmPtr warm;           // for the children (null without the active-set method)
  std::vector<double> lb, ub;  // the bounds it was solved with (the safe bound and d refer to them)
};

struct Node {
  double bound = -kInf;  // lower bound on the objective in the node (from its parent)
  int depth = 0;
  long long id = 0;
  std::vector<BoundChange> changes;  // from the root bounds, applied in order
  // Pseudocost bookkeeping: the branching that created the node.
  Index var = -1;
  int dir = 0;              // -1 down, +1 up
  double dist = 0.0;        // distance of the parent's value to the new bound
  double parent_obj = 0.0;  // the parent's relaxation objective
  WarmPtr warm;             // the parent's working set (shared by both children)
};

struct NodeOrder {
  bool operator()(const Node& a, const Node& b) const {  // true: a after b (min-heap on bound)
    if (a.bound != b.bound) return a.bound > b.bound;
    if (a.depth != b.depth) return a.depth < b.depth;
    return a.id > b.id;
  }
};

// Memory of a queued node (a shared warm start is counted half per child).
long long node_bytes(const Node& n) {
  return static_cast<long long>(sizeof(Node) + n.changes.size() * sizeof(BoundChange)) +
         (n.warm ? n.warm->bytes() / 2 : 0);
}

double frac(double v) { return std::abs(v - std::round(v)); }

class Miqp {
 public:
  Miqp(const MipProblem& p, const SolverParams& params, const Logger& log, const Timer& timer, const MiqpOptions& opt)
      : p_(p), params_(params), log_(log), timer_(timer), opt_(opt), prop_(p) {}

  MiqpResult run();

 private:
  QpSolution relax(const std::vector<double>& lb, const std::vector<double>& ub, bool root = false);
  // The node relaxation: the active-set method from `warm` when available, else (or if it
  // fails) the barrier.
  NodeQp solve_node(const std::vector<double>& lb, const std::vector<double>& ub, const WarmPtr& warm,
                    bool fallback = true);
  NodeQp from_barrier(const QpSolution& s, const std::vector<double>& lb, const std::vector<double>& ub);
  // Reduced-cost tightening of the integer columns' bounds from a safe bound and its reduced
  // costs: appends the changes, returns their number.
  // slb / sub: the bounds the safe bound was computed with; lb / ub: the current ones (at least
  // as tight). Only changes that tighten lb / ub are appended.
  int rc_tighten(double bound, const std::vector<double>& d, const std::vector<double>& slb,
                 const std::vector<double>& sub, const std::vector<double>& lb, const std::vector<double>& ub,
                 std::vector<BoundChange>& out) const;
  bool propagate(std::vector<double>& lb, std::vector<double>& ub, const std::vector<Index>& seeds);
  // Integer columns of x rounded; accepted (true) if feasible, else the continuous columns are
  // re-optimized with the integers fixed. Updates the incumbent if better.
  bool try_point(const std::vector<double>& x, const char* source, const WarmPtr& warm = nullptr);
  bool improves(double obj) const {
    return !is_finite(inc_obj_) || obj < inc_obj_ - 1e-9 * std::max(1.0, std::abs(inc_obj_));
  }
  double cutoff() const;  // nodes with a bound >= this are pruned
  void qp_dive(std::vector<double> lb, std::vector<double> ub, std::vector<double> x, WarmPtr warm, int max_solves);
  void run_fj(const std::vector<double>& xroot);  // xroot empty: feasibility only
  void simple_rounding(const std::vector<double>& x);
  bool limits_hit();
  double global_bound(double current) const;
  void log_line(double current_bound, bool force);
  double display(double v) const { return p_.display_sense * v; }

  const MipProblem& p_;
  const SolverParams& params_;
  const Logger& log_;
  const Timer& timer_;
  const MiqpOptions& opt_;
  DomainPropagator prop_;
  bool has_cont_ = false;
  std::unique_ptr<ActiveSetQp> as_;
  long long as_budget_ = 0;
  long long as_solves_ = 0, as_fail_ = 0, as_iter_ = 0, rc_fixed_ = 0, as_flops_ = 0;
  bool sb_enabled_ = true;  // reliability branching (with the active-set method)
  long long sb_solves_ = 0, sb_fixed_ = 0;

  std::vector<double> inc_x_;
  double inc_obj_ = kInf;
  std::vector<double> pc_sum_[2];
  std::vector<int> pc_cnt_[2];
  std::priority_queue<Node, std::vector<Node>, NodeOrder> open_;
  long long nodes_ = 0, next_id_ = 0, barrier_it_ = 0;
  double pruned_bound_ = kInf;  // smallest bound of a node pruned by the gap (below the incumbent)
  double lost_bound_ = kInf;    // smallest bound of a node whose relaxation could not be solved
  Status limit_status_ = Status::Unknown;
  double last_log_ = -1.0;
  long long open_bytes_ = 0;
};

QpSolution Miqp::relax(const std::vector<double>& lb, const std::vector<double>& ub, bool root) {
  QpOptions o;
  o.time_limit = params_.time_limit - timer_.elapsed();
  o.threads = root || p_.num_cols() > 5000 ? opt_.threads : 1;
  o.log = root && log_.verbosity() >= 2 ? &log_ : nullptr;
  QpSolution s = solve_qp(p_, lb, ub, o, opt_.diag_shift);
  barrier_it_ += s.iterations;
  return s;
}

NodeQp Miqp::from_barrier(const QpSolution& s, const std::vector<double>& lb, const std::vector<double>& ub) {
  NodeQp r;
  r.status = s.status;
  if (s.status != QpStatus::Optimal) return r;
  r.x = s.x;
  r.objective = s.objective;
  r.bound = s.dual_bound;
  if (as_) {
    std::vector<double> qx;
    as_->q_times(s.x, qx);
    SafeBound sb = qp_safe_bound(p_, qx, s.x, s.y, lb, ub);
    if (is_finite(sb.bound)) {
      r.bound = sb.bound;
      r.d = std::move(sb.d);
    }
    auto w = std::make_shared<Warm>();
    w->ws = as_->guess(s.x, lb, ub);
    w->x = s.x;
    r.warm = std::move(w);
  }
  return r;
}

NodeQp Miqp::solve_node(const std::vector<double>& lb, const std::vector<double>& ub, const WarmPtr& warm,
                        bool fallback) {
  if (as_ && warm) {
    ++as_solves_;
    ActiveSetResult a = as_->solve(lb, ub, warm->ws, warm->x, as_budget_, warm->state);
    as_iter_ += a.iterations;
    as_flops_ += a.flops;
    if (a.status == ActiveSetResult::Status::Infeasible) {
      NodeQp r;
      r.status = QpStatus::Infeasible;
      return r;
    }
    if (a.status == ActiveSetResult::Status::Optimal && is_finite(a.safe.bound)) {
      NodeQp r;
      r.status = QpStatus::Optimal;
      r.objective = a.objective;
      r.bound = a.safe.bound;
      r.d = std::move(a.safe.d);
      auto w = std::make_shared<Warm>();
      w->ws = std::move(a.ws);
      w->state = a.state;
      w->x = a.x;
      r.x = std::move(a.x);
      r.warm = std::move(w);
      r.lb = lb;
      r.ub = ub;
      return r;
    }
    if (!fallback) {
      NodeQp r;
      r.lb = lb;
      r.ub = ub;
      return r;
    }
    ++as_fail_;
    log_.log(3, "MIQP:   active-set solve failed (%s), barrier\n",
             a.status == ActiveSetResult::Status::Optimal ? "no safe bound" : a.reason);
  }
  NodeQp r = from_barrier(relax(lb, ub), lb, ub);
  r.lb = lb;
  r.ub = ub;
  return r;
}

int Miqp::rc_tighten(double bound, const std::vector<double>& d, const std::vector<double>& slb,
                     const std::vector<double>& sub, const std::vector<double>& lb, const std::vector<double>& ub,
                     std::vector<BoundChange>& out) const {
  // f(x) >= bound + d_j (x_j - slb_j) for d_j > 0 (and the mirror image for d_j < 0): values of x_j
  // beyond room / |d_j| from that bound cannot improve on the incumbent.
  if (d.empty() || slb.empty() || !is_finite(inc_obj_) || !is_finite(bound)) return 0;
  const double room = cutoff() - bound;
  if (!(room >= 0.0)) return 0;
  int cnt = 0;
  for (Index j = 0; j < p_.num_cols(); ++j) {
    if (!p_.is_int[j] || lb[j] == ub[j]) continue;
    const double dj = d[j];
    if (dj > 0.0 && is_finite(slb[j])) {
      const double nu = slb[j] + std::floor(room / dj + 1e-6);
      if (nu < ub[j]) out.push_back({j, lb[j], nu}), ++cnt;
    } else if (dj < 0.0 && is_finite(sub[j])) {
      const double nl = sub[j] - std::floor(room / -dj + 1e-6);
      if (nl > lb[j]) out.push_back({j, nl, ub[j]}), ++cnt;
    }
  }
  return cnt;
}

bool Miqp::propagate(std::vector<double>& lb, std::vector<double>& ub, const std::vector<Index>& seeds) {
  std::vector<Index> tightened;
  return prop_.propagate(lb, ub, seeds, tightened);
}

double Miqp::cutoff() const {
  if (!is_finite(inc_obj_)) return kInf;
  const double tol = std::max(params_.mip_abs_gap, params_.mip_rel_gap * std::abs(inc_obj_));
  return inc_obj_ - tol;
}

bool Miqp::try_point(const std::vector<double>& x_in, const char* source, const WarmPtr& warm) {
  const Index n = p_.num_cols();
  std::vector<double> x = x_in;
  for (Index j = 0; j < n; ++j) {
    if (p_.is_int[j]) x[j] = std::round(x[j]);
    x[j] = std::min(std::max(x[j], p_.col_lb[j]), p_.col_ub[j]);
  }
  double int_viol = 0.0;
  if (!(p_.max_violation(x, &int_viol) <= kFeasTol)) {
    log_.log(3, "MIQP:   rounded point (%s) violates by %.3g\n", source, p_.max_violation(x));
    if (!has_cont_) return false;
    std::vector<double> lb = p_.col_lb, ub = p_.col_ub;
    for (Index j = 0; j < n; ++j)
      if (p_.is_int[j]) lb[j] = ub[j] = x[j];
    const NodeQp s = solve_node(lb, ub, warm);
    log_.log(3, "MIQP:   completion QP: %s\n", to_string(s.status));
    if (s.status != QpStatus::Optimal) return false;
    x = s.x;
    for (Index j = 0; j < n; ++j)
      if (p_.is_int[j]) x[j] = lb[j];
    if (!(p_.max_violation(x) <= kFeasTol)) return false;
  }
  const double obj = p_.objective(x);
  log_.log(3, "MIQP:   candidate (%s) objective %.10g\n", source, display(obj));
  if (!improves(obj)) return true;
  inc_obj_ = obj;
  inc_x_ = std::move(x);
  log_.log(2, "MIQP: new incumbent %.10g (%s) at node %lld, %.2fs\n", display(obj), source, nodes_, timer_.elapsed());
  return true;
}

void Miqp::simple_rounding(const std::vector<double>& x) {
  // Pure-integer problems: rounding is a complete candidate (no QP needed).
  std::vector<double> r = x;
  for (Index j = 0; j < p_.num_cols(); ++j)
    if (p_.is_int[j]) r[j] = std::round(r[j]);
  if (p_.max_violation(r) <= kFeasTol && improves(p_.objective(r))) try_point(r, "rounding");
}

void Miqp::qp_dive(std::vector<double> lb, std::vector<double> ub, std::vector<double> x, WarmPtr warm,
                   int max_solves) {
  const Index n = p_.num_cols();
  int solves = 0;
  while (solves < max_solves && !limits_hit()) {
    // Least fractional unfixed integer column; integral ones are fixed along with it.
    Index best = -1;
    double best_f = kInf;
    for (Index j = 0; j < n; ++j) {
      if (!p_.is_int[j] || lb[j] == ub[j]) continue;
      const double f = frac(x[j]);
      if (f > kIntTol && f < best_f) best_f = f, best = j;
    }
    if (best < 0) {
      try_point(x, "QP dive", warm);
      return;
    }
    bool ok = false;
    for (int attempt = 0; attempt < 2 && !ok; ++attempt) {
      // Attempt 0: the rounded value plus the integral columns; attempt 1 (backtrack): only the
      // other direction.
      double v = std::round(x[best]);
      if (attempt == 1) v = v > x[best] ? v - 1.0 : v + 1.0;
      if (v < lb[best] || v > ub[best]) continue;
      std::vector<double> l2 = lb, u2 = ub;
      std::vector<Index> seeds;
      for (Index j = 0; j < n && attempt == 0; ++j) {
        if (!p_.is_int[j] || l2[j] == u2[j] || j == best || frac(x[j]) > kIntTol) continue;
        l2[j] = u2[j] = std::round(x[j]);
        seeds.push_back(j);
      }
      l2[best] = u2[best] = v;
      seeds.push_back(best);
      if (!propagate(l2, u2, seeds)) continue;
      NodeQp s = solve_node(l2, u2, warm);
      ++solves;
      if (s.status != QpStatus::Optimal) continue;
      if (s.bound >= cutoff()) return;
      lb = std::move(l2), ub = std::move(u2), x = std::move(s.x);
      if (s.warm && warm) warm = std::move(s.warm);
      ok = true;
    }
    if (!ok) return;
  }
}

void Miqp::run_fj(const std::vector<double>& xroot) {
  // The MILP with the objective linearized at the root relaxation's point (gradient c + Q x).
  MipProblem lin = p_;
  lin.Q = SparseMatrix();
  if (xroot.empty()) std::fill(lin.cost.begin(), lin.cost.end(), 0.0);  // feasibility only
  for (Index j = 0; j < p_.Q.num_major && !xroot.empty(); ++j) {
    for (Index k = p_.Q.begin(j); k < p_.Q.end(j); ++k) {
      const Index i = p_.Q.index[k];
      lin.cost[i] += p_.Q.value[k] * xroot[j];
      if (i != j) lin.cost[j] += p_.Q.value[k] * xroot[i];
    }
  }
  FjLimits lim;
  lim.effort = std::min<long long>(20000000, 200LL * (p_.A.nnz() + p_.num_cols()) + 100000);
  lim.stall_effort = lim.effort / 2;
  lim.time_limit = std::max(0.0, std::min(10.0, 0.1 * (params_.time_limit - timer_.elapsed())));
  lim.seed = params_.seed;
  std::vector<std::vector<double>> found;
  feasibility_jump(lin, lim, [&](const std::vector<double>& x, double) {
    if (found.size() >= 4) found.erase(found.begin());
    found.push_back(x);
  });
  for (auto it = found.rbegin(); it != found.rend(); ++it) try_point(*it, "feasibility jump");
}

bool Miqp::limits_hit() {
  if (limit_status_ != Status::Unknown) return true;
  if (timer_.elapsed() >= params_.time_limit) limit_status_ = Status::TimeLimit;
  else if (params_.node_limit >= 0 && nodes_ >= params_.node_limit) limit_status_ = Status::NodeLimit;
  else if (params_.mem_limit_mb > 0 && open_bytes_ > params_.mem_limit_mb * 1024LL * 1024LL) limit_status_ = Status::MemoryLimit;
  return limit_status_ != Status::Unknown;
}

double Miqp::global_bound(double current) const {
  double b = std::min({current, pruned_bound_, lost_bound_});
  if (!open_.empty()) b = std::min(b, open_.top().bound);
  return std::min(b, inc_obj_);
}

void Miqp::log_line(double current_bound, bool force) {
  const double t = timer_.elapsed();
  if (!force && t - last_log_ < 5.0) return;
  last_log_ = t;
  const double gb = global_bound(current_bound);
  double gap = kInf;
  if (is_finite(inc_obj_) && is_finite(gb)) gap = std::abs(inc_obj_ - gb) / std::max(std::abs(inc_obj_), 1e-10);
  log_.log(1, "MIQP: %8lld nodes, %7zu open, incumbent %14.8g, bound %14.8g, gap %7.3f%%, %.1fs\n", nodes_,
           open_.size(), is_finite(inc_obj_) ? display(inc_obj_) : display(kInf), display(gb), 100.0 * gap, t);
}

MiqpResult Miqp::run() {
  MiqpResult res;
  const Index n = p_.num_cols();
  for (Index j = 0; j < n; ++j) has_cont_ = has_cont_ || !p_.is_int[j];
  for (int d = 0; d < 2; ++d) {
    pc_sum_[d].assign(static_cast<size_t>(n), 0.0);
    pc_cnt_[d].assign(static_cast<size_t>(n), 0);
  }
  Index nint = 0;
  for (Index j = 0; j < n; ++j) nint += p_.is_int[j] ? 1 : 0;
  log_.log(1, "MIQP: %d rows, %d columns (%d integer), branch-and-bound on QP relaxations\n", p_.num_rows(), n, nint);

  // Root.
  std::vector<double> rlb = p_.col_lb, rub = p_.col_ub;
  std::vector<Index> all(static_cast<size_t>(n));
  for (Index j = 0; j < n; ++j) all[j] = j;
  auto finish = [&](Status st) {
    res.status = st;
    res.nodes = nodes_;
    res.barrier_iterations = barrier_it_;
    res.lp_iterations = as_iter_;
    res.has_solution = is_finite(inc_obj_);
    if (res.has_solution) {
      res.x = inc_x_;
      res.objective = inc_obj_;
    }
    return res;
  };
  if (!propagate(rlb, rub, all)) {
    res.bound = kInf;
    return finish(Status::Infeasible);
  }
  if (!opt_.start.empty() && static_cast<Index>(opt_.start.size()) == n) {
    bool complete = true, ints = true;
    for (Index j = 0; j < n; ++j) {
      if (std::isnan(opt_.start[j])) {
        complete = false;
        if (p_.is_int[j]) ints = false;
      }
    }
    if (complete || ints) {
      std::vector<double> s = opt_.start;
      for (Index j = 0; j < n; ++j)
        if (std::isnan(s[j])) s[j] = std::max(p_.col_lb[j], std::min(0.0, p_.col_ub[j]));
      const bool ok = try_point(s, "start");
      log_.log(1, "MIQP: start solution %s\n", ok ? "accepted" : "infeasible, ignored");
    } else {
      log_.log(1, "MIQP: partial start solutions without all integer values are ignored\n");
    }
  }
  const QpSolution root = relax(rlb, rub, true);
  ++nodes_;
  if (root.status == QpStatus::Infeasible) {
    res.bound = kInf;
    return finish(Status::Infeasible);
  }
  if (root.status == QpStatus::Unbounded || root.status == QpStatus::InfeasibleOrUnbounded) {
    // An unbounded relaxation of a MIQP with linear constraints: unbounded if feasible (a rational
    // recession direction with Q d = 0, c^T d < 0 scales to an integral one). Look for a point.
    log_.log(1, "MIQP: the root relaxation is unbounded\n");
    if (root.status == QpStatus::Unbounded && !is_finite(inc_obj_)) run_fj({});
    return finish(is_finite(inc_obj_) ? Status::Unbounded : Status::InfeasibleOrUnbounded);
  }
  if (root.status == QpStatus::TimeLimit) {
    res.bound = -kInf;
    return finish(Status::TimeLimit);
  }
  if (root.status != QpStatus::Optimal) {
    log_.log(1, "MIQP: the root relaxation failed (%s)\n", to_string(root.status));
    return finish(Status::NumericalError);
  }
  log_.log(1, "MIQP: root relaxation %.10g (%d barrier iterations, %.2fs)\n", display(root.objective), root.iterations,
           timer_.elapsed());
  // Node relaxations by the dual active-set method (dense; warm started from the parent) when
  // the root's working set can be completed by it within its work limit.
  NodeQp rootq;
  if (n <= 3000 && static_cast<double>(n) * n * 8.0 <= 0.25 * params_.mem_limit_mb * 1048576.0) {
    as_ = std::make_unique<ActiveSetQp>(p_, opt_.diag_shift);
    const double nn = static_cast<double>(n);
    as_budget_ = static_cast<long long>(std::min(4e9, 5.0 * nn * nn * nn + 2e6));
    auto w0 = std::make_shared<Warm>();
    w0->ws = as_->guess(root.x, rlb, rub);
    w0->x = root.x;
    const long long budget = as_budget_;
    as_budget_ = 4 * budget;
    rootq = solve_node(rlb, rub, w0);
    as_budget_ = budget;
    if (as_fail_ > 0 || rootq.status != QpStatus::Optimal) {
      log_.log(1, "MIQP: the active-set method failed at the root; barrier relaxations\n");
      as_.reset();
      rootq = NodeQp();
    } else {
      log_.log(1, "MIQP: root active set: %zu entries, %lld iterations, bound %.10g\n", rootq.warm->ws.size(), as_iter_,
               display(rootq.bound));
    }
  }
  if (rootq.status != QpStatus::Optimal) {
    rootq.status = QpStatus::Optimal;
    rootq.x = root.x;
    rootq.objective = root.objective;
    rootq.bound = root.dual_bound;
  }
  res.root_bound = rootq.bound;
  // Root heuristics.
  // (on the barrier's central point: rounding it works better than rounding a vertex-like one)
  if (!has_cont_) simple_rounding(root.x);
  if (has_cont_) try_point(root.x, "rounding", rootq.warm);
  if (!limits_hit()) run_fj(root.x);
  if (!limits_hit()) qp_dive(rlb, rub, root.x, nullptr, std::min<int>(200, 2 * nint + 2));
  if (is_finite(inc_obj_)) {
    log_.log(1, "MIQP: root heuristics: incumbent %.10g (%.2fs)\n", display(inc_obj_), timer_.elapsed());
  }
  // Root reduced-cost fixing (global), redone whenever the incumbent improves.
  const std::vector<double> root_d = rootq.d;
  const double root_lb = rootq.bound;
  const std::vector<double> root_slb = rootq.lb, root_sub = rootq.ub;
  double root_fix_inc = kInf;
  auto root_fixing = [&]() {
    if (root_d.empty() || !(inc_obj_ < root_fix_inc)) return;
    root_fix_inc = inc_obj_;
    std::vector<BoundChange> ch;
    if (rc_tighten(root_lb, root_d, root_slb, root_sub, rlb, rub, ch) == 0) return;
    std::vector<Index> seeds;
    for (const BoundChange& c : ch) {
      rlb[c.col] = std::max(rlb[c.col], c.lb);
      rub[c.col] = std::min(rub[c.col], c.ub);
      seeds.push_back(c.col);
    }
    rc_fixed_ += static_cast<long long>(ch.size());
    log_.log(2, "MIQP: root reduced-cost tightening: %zu columns\n", ch.size());
    // A failed propagation means no improving point exists; the search then prunes everything.
    bool empty = false;
    for (const BoundChange& c : ch) empty = empty || rlb[c.col] > rub[c.col];
    if (empty || !propagate(rlb, rub, seeds)) {
      for (Index j = 0; j < n; ++j) rub[j] = std::min(rub[j], rlb[j] - 1.0);
    }
  };
  root_fixing();

  // Branch-and-bound. The current node of a plunge is held outside the queue.
  Node cur;
  cur.bound = rootq.bound;
  cur.warm = rootq.warm;
  bool have_cur = true;
  NodeQp cur_sol = rootq;
  bool cur_solved = true;  // the root's relaxation is already solved
  long long next_dive = 50;
  while (true) {
    if (!have_cur) {
      if (open_.empty()) break;
      cur = open_.top();
      open_.pop();
      open_bytes_ -= node_bytes(cur);
      have_cur = true;
      cur_solved = false;
    }
    root_fixing();
    if (cur.bound >= cutoff()) {
      if (cur.bound < inc_obj_) pruned_bound_ = std::min(pruned_bound_, cur.bound);
      have_cur = false;
      continue;
    }
    {
      // Gap reached?
      const double gb = global_bound(cur.bound);
      if (is_finite(inc_obj_) &&
          (inc_obj_ - gb <= params_.mip_abs_gap || inc_obj_ - gb <= params_.mip_rel_gap * std::abs(inc_obj_))) {
        pruned_bound_ = std::min(pruned_bound_, cur.bound);
        break;
      }
    }
    if (limits_hit()) break;
    log_line(cur.bound, false);
    std::vector<double> lb = rlb, ub = rub;
    std::vector<Index> seeds;
    bool infeasible = false;
    for (const BoundChange& c : cur.changes) {
      lb[c.col] = std::max(lb[c.col], c.lb);
      ub[c.col] = std::min(ub[c.col], c.ub);
      seeds.push_back(c.col);
    }
    for (Index j = 0; j < n && !infeasible; ++j) infeasible = lb[j] > ub[j];
    if (!cur_solved) {
      ++nodes_;
      if (infeasible || !propagate(lb, ub, seeds)) {
        have_cur = false;
        continue;
      }
      cur_sol = solve_node(lb, ub, cur.warm);
      if (cur_sol.status == QpStatus::Infeasible) {
        have_cur = false;
        continue;
      }
      if (cur_sol.status == QpStatus::TimeLimit) {
        limit_status_ = Status::TimeLimit;
        break;
      }
    } else if (infeasible) {
      have_cur = false;
      continue;
    }
    cur_solved = false;
    if (cur_sol.status == QpStatus::Infeasible) {
      have_cur = false;
      continue;
    }
    if (cur_sol.status != QpStatus::Optimal) {
      // Relaxation not solved: split the domain of an unfixed integer column (no bound).
      Index j = -1;
      for (Index k = 0; k < n && j < 0; ++k)
        if (p_.is_int[k] && lb[k] < ub[k]) j = k;
      if (j < 0) {
        lost_bound_ = std::min(lost_bound_, cur.bound);
        have_cur = false;
        continue;
      }
      const double mid = is_finite(lb[j]) && is_finite(ub[j]) ? std::floor(0.5 * (lb[j] + ub[j]))
                         : is_finite(lb[j])                    ? lb[j]
                         : is_finite(ub[j])                    ? ub[j] - 1.0
                                                               : 0.0;
      for (int side = 0; side < 2; ++side) {
        Node c;
        c.bound = cur.bound;
        c.depth = cur.depth + 1;
        c.id = next_id_++;
        c.changes = cur.changes;
        c.changes.push_back(side == 0 ? BoundChange{j, -kInf, mid} : BoundChange{j, mid + 1.0, kInf});
        c.warm = cur.warm;
        open_bytes_ += node_bytes(c);
        open_.push(std::move(c));
      }
      have_cur = false;
      continue;
    }
    // The box the relaxation was solved with (propagated), within the current one.
    if (!cur_sol.lb.empty()) {
      for (Index j = 0; j < n; ++j) {
        lb[j] = std::max(lb[j], cur_sol.lb[j]);
        ub[j] = std::min(ub[j], cur_sol.ub[j]);
      }
    }
    // Reduced-cost tightening: re-solve the node with the tighter bounds (a few rounds).
    for (int round = 0; round < 3 && cur_sol.status == QpStatus::Optimal; ++round) {
      const double nb = std::max(cur.bound, cur_sol.bound);
      if (nb >= cutoff()) break;
      std::vector<BoundChange> ch;
      // The safe bound itself (not the parent's) goes with its reduced costs and box.
      if (rc_tighten(cur_sol.bound, cur_sol.d, cur_sol.lb, cur_sol.ub, lb, ub, ch) == 0) break;
      rc_fixed_ += static_cast<long long>(ch.size());
      std::vector<Index> s2;
      bool empty = false;
      for (const BoundChange& c : ch) {
        lb[c.col] = std::max(lb[c.col], c.lb);
        ub[c.col] = std::min(ub[c.col], c.ub);
        empty = empty || lb[c.col] > ub[c.col];
        s2.push_back(c.col);
        cur.changes.push_back(c);
      }
      if (empty || !propagate(lb, ub, s2)) {
        cur_sol.status = QpStatus::Infeasible;
        break;
      }
      const WarmPtr w = cur_sol.warm ? cur_sol.warm : cur.warm;
      cur.bound = nb;  // still valid for the tightened node
      cur_sol = solve_node(lb, ub, w);
    }
    if (cur_sol.status == QpStatus::Infeasible) {
      have_cur = false;
      continue;
    }
    if (cur_sol.status == QpStatus::TimeLimit) {
      limit_status_ = Status::TimeLimit;
      break;
    }
    if (cur_sol.status != QpStatus::Optimal) {
      // The re-solve failed: keep the node's bound and split it (as above, via the queue).
      Node c = cur;
      c.id = next_id_++;
      Index j = -1;
      for (Index k = 0; k < n && j < 0; ++k)
        if (p_.is_int[k] && lb[k] < ub[k]) j = k;
      if (j < 0) {
        lost_bound_ = std::min(lost_bound_, cur.bound);
        have_cur = false;
        continue;
      }
      const double mid = is_finite(lb[j]) && is_finite(ub[j]) ? std::floor(0.5 * (lb[j] + ub[j]))
                         : is_finite(lb[j])                    ? lb[j]
                         : is_finite(ub[j])                    ? ub[j] - 1.0
                                                               : 0.0;
      for (int side = 0; side < 2; ++side) {
        Node s = c;
        s.id = next_id_++;
        s.depth = cur.depth + 1;
        s.changes.push_back(side == 0 ? BoundChange{j, -kInf, mid} : BoundChange{j, mid + 1.0, kInf});
        open_bytes_ += node_bytes(s);
        open_.push(std::move(s));
      }
      have_cur = false;
      continue;
    }
    const double node_obj = cur_sol.objective;
    const double node_bound = std::max(cur.bound, cur_sol.bound);
    log_.log(3, "MIQP: node %lld depth %d: relaxation %.10g (bound %.10g), %zu bound changes\n", nodes_, cur.depth,
             display(node_obj), display(node_bound), cur.changes.size());
    // Pseudocost update (once per node).
    if (cur.var >= 0 && cur.dist > 0.0) {
      const int d = cur.dir > 0 ? 1 : 0;
      pc_sum_[d][cur.var] += std::max(0.0, node_obj - cur.parent_obj) / cur.dist;
      ++pc_cnt_[d][cur.var];
      cur.var = -1;
    }
    if (node_bound >= cutoff()) {
      if (node_bound < inc_obj_) pruned_bound_ = std::min(pruned_bound_, node_bound);
      have_cur = false;
      continue;
    }
    const std::vector<double>& x = cur_sol.x;
    // Branching candidates: pseudocost scores (product of the estimated gains).
    Index best = -1;
    double best_score = -1.0;
    std::vector<std::pair<double, Index>> cands;
    {
      double avg[2] = {0.0, 0.0};
      int cnt[2] = {0, 0};
      for (Index j = 0; j < n; ++j) {
        for (int d = 0; d < 2; ++d) {
          if (pc_cnt_[d][j] > 0) avg[d] += pc_sum_[d][j] / pc_cnt_[d][j], ++cnt[d];
        }
      }
      const bool use_pc = cnt[0] > 0 && cnt[1] > 0;
      for (int d = 0; d < 2; ++d) avg[d] = cnt[d] > 0 ? avg[d] / cnt[d] : 1.0;
      for (Index j = 0; j < n; ++j) {
        if (!p_.is_int[j] || lb[j] == ub[j]) continue;
        const double f = x[j] - std::floor(x[j]);
        if (frac(x[j]) <= kIntTol) continue;
        double score;
        if (use_pc) {
          const double pd = pc_cnt_[0][j] > 0 ? pc_sum_[0][j] / pc_cnt_[0][j] : avg[0];
          const double pu = pc_cnt_[1][j] > 0 ? pc_sum_[1][j] / pc_cnt_[1][j] : avg[1];
          score = std::max(f * pd, 1e-6) * std::max((1.0 - f) * pu, 1e-6);
        } else {
          score = std::min(f, 1.0 - f);
        }
        cands.push_back({score, j});
        if (score > best_score + 1e-12) best_score = score, best = j;
      }
    }
    // Reliability branching: strong branching (active-set solves from this node) on the best
    // candidates whose pseudocosts rest on fewer than kRel observations. A child that is
    // infeasible or cut off fixes the column to the other side (the node continues as that
    // child, already solved).
    if (as_ && cur_sol.warm && best >= 0 && sb_enabled_) {
      constexpr int kRel = 4, kMaxEval = 8, kLookahead = 4;
      std::stable_sort(cands.begin(), cands.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
      double sb_best = -1.0;
      Index sb_j = -1;
      int evaluated = 0, since = 0;
      bool fixed = false, node_done = false;
      for (const auto& [pscore, j] : cands) {
        if (evaluated >= kMaxEval || since >= kLookahead || limits_hit()) break;
        if (std::min(pc_cnt_[0][j], pc_cnt_[1][j]) >= kRel) {
          if (pscore > sb_best) sb_best = pscore, sb_j = j, since = 0;
          else ++since;
          continue;
        }
        ++evaluated;
        const double v = x[j];
        NodeQp child[2];
        double gain[2];
        bool cut[2];
        std::vector<double> box_lb[2], box_ub[2];
        for (int dir = 0; dir < 2; ++dir) {
          std::vector<double>& l2 = box_lb[dir];
          std::vector<double>& u2 = box_ub[dir];
          l2 = lb, u2 = ub;
          if (dir == 0) u2[j] = std::floor(v);
          else l2[j] = std::ceil(v);
          const std::vector<Index> sd{j};
          ++sb_solves_;
          if (l2[j] > u2[j] || !propagate(l2, u2, sd)) {
            child[dir].status = QpStatus::Infeasible;
          } else {
            child[dir] = solve_node(l2, u2, cur_sol.warm, /*fallback=*/false);
          }
          cut[dir] = child[dir].status == QpStatus::Infeasible ||
                     (child[dir].status == QpStatus::Optimal && child[dir].bound >= cutoff());
          gain[dir] = child[dir].status == QpStatus::Optimal ? std::max(0.0, child[dir].objective - node_obj) : -1.0;
          if (gain[dir] >= 0.0) {
            const double dist = dir == 0 ? v - std::floor(v) : std::ceil(v) - v;
            pc_sum_[dir][j] += gain[dir] / dist;
            ++pc_cnt_[dir][j];
          }
        }
        if (cut[0] && cut[1]) {
          for (int dir = 0; dir < 2; ++dir)
            if (child[dir].status == QpStatus::Optimal && child[dir].bound < inc_obj_)
              pruned_bound_ = std::min(pruned_bound_, child[dir].bound);
          node_done = true;
          break;
        }
        if (cut[0] || cut[1]) {
          const int keep = cut[0] ? 1 : 0;
          if (cut[1 - keep] && child[1 - keep].status == QpStatus::Optimal && child[1 - keep].bound < inc_obj_)
            pruned_bound_ = std::min(pruned_bound_, child[1 - keep].bound);
          cur.changes.push_back(keep == 0 ? BoundChange{j, -kInf, std::floor(v)} : BoundChange{j, std::ceil(v), kInf});
          ++sb_fixed_;
          if (child[keep].status == QpStatus::Optimal) {
            cur.bound = node_bound;
            cur_sol = std::move(child[keep]);
            fixed = true;
          } else {
            // The remaining child's active-set solve failed: solve it as a node (barrier fallback).
            cur.bound = node_bound;
            const WarmPtr w = cur_sol.warm;
            cur_sol = solve_node(box_lb[keep], box_ub[keep], w);
            fixed = true;
          }
          break;
        }
        const double sc = std::max(gain[0], 1e-6) * std::max(gain[1], 1e-6);
        if (gain[0] >= 0.0 && gain[1] >= 0.0 && sc > sb_best) sb_best = sc, sb_j = j, since = 0;
        else ++since;
      }
      if (node_done) {
        have_cur = false;
        continue;
      }
      if (fixed) {
        cur_solved = true;  // continue with the remaining child as this node
        continue;
      }
      if (sb_j >= 0) best = sb_j;
    }
    if (log_.verbosity() >= 3) {
      std::string xs;
      for (Index j = 0; j < n; ++j) xs += " " + std::to_string(x[j]) + (p_.is_int[j] ? "i" : "");
      log_.log(3, "MIQP:   x =%s; branching column %d\n", xs.c_str(), best);
    }
    const WarmPtr warm = cur_sol.warm ? cur_sol.warm : cur.warm;
    if (best < 0) {
      // Integral relaxation: a candidate solution, and the node is done if it is feasible.
      if (try_point(x, "relaxation", warm)) {
        have_cur = false;
        continue;
      }
      // The rounded point cannot be completed: split three ways on the column farthest from
      // integrality (or any unfixed integer column).
      Index j = -1;
      double fj = -1.0;
      for (Index k = 0; k < n; ++k) {
        if (p_.is_int[k] && lb[k] < ub[k] && frac(x[k]) > fj) fj = frac(x[k]), j = k;
      }
      if (j < 0) {
        have_cur = false;  // all integers fixed: the node's QP is exactly the completion
        continue;
      }
      const double v = std::round(x[j]);
      const BoundChange parts[3] = {{j, -kInf, v - 1.0}, {j, v, v}, {j, v + 1.0, kInf}};
      for (const BoundChange& bc : parts) {
        if (std::max(bc.lb, lb[j]) > std::min(bc.ub, ub[j])) continue;
        Node c;
        c.bound = node_bound;
        c.depth = cur.depth + 1;
        c.id = next_id_++;
        c.changes = cur.changes;
        c.changes.push_back(bc);
        c.warm = warm;
        open_bytes_ += node_bytes(c);
        open_.push(std::move(c));
      }
      have_cur = false;
      continue;
    }
    if (!has_cont_) simple_rounding(x);
    if (nodes_ >= next_dive && !limits_hit()) {
      next_dive = 2 * nodes_ + 50;
      qp_dive(lb, ub, x, warm, std::min<int>(100, nint + 2));
      if (node_bound >= cutoff()) {
        have_cur = false;
        continue;
      }
    }
    // Branch: both children; continue with the one in the rounding direction.
    const double v = x[best];
    Node down, up;
    for (Node* c : {&down, &up}) {
      c->bound = node_bound;
      c->depth = cur.depth + 1;
      c->id = next_id_++;
      c->changes = cur.changes;
      c->var = best;
      c->parent_obj = node_obj;
      c->warm = warm;
    }
    down.changes.push_back({best, -kInf, std::floor(v)});
    down.dir = -1;
    down.dist = v - std::floor(v);
    up.changes.push_back({best, std::ceil(v), kInf});
    up.dir = 1;
    up.dist = std::ceil(v) - v;
    // Plunge into the child with the smaller estimated degradation (pseudocosts) once both
    // directions have observations, else the rounding direction.
    bool go_up = v - std::floor(v) >= 0.5;
    if (pc_cnt_[0][best] > 0 && pc_cnt_[1][best] > 0) {
      const double ed = (v - std::floor(v)) * pc_sum_[0][best] / pc_cnt_[0][best];
      const double eu = (std::ceil(v) - v) * pc_sum_[1][best] / pc_cnt_[1][best];
      go_up = eu < ed;
    }
    Node& next = go_up ? up : down;
    Node& other = go_up ? down : up;
    open_bytes_ += node_bytes(other);
    open_.push(std::move(other));
    cur = std::move(next);
    have_cur = true;
    cur_solved = false;
  }
  const double cur_bound = have_cur ? cur.bound : kInf;
  double bound = global_bound(cur_bound);
  log_line(cur_bound, true);
  if (as_) {
    log_.log(1, "MIQP: active-set solves %lld (%lld iterations, %.3g flops, %lld fell back to the barrier), %lld reduced-cost "
             "tightenings, %lld strong-branching solves (%lld fixings)\n", as_solves_, as_iter_,
             static_cast<double>(as_flops_), as_fail_, rc_fixed_, sb_solves_, sb_fixed_);
  }
  if (limit_status_ != Status::Unknown) {
    res.bound = bound;
    return finish(limit_status_);
  }
  // Search complete (or gap reached).
  if (!is_finite(inc_obj_)) {
    if (is_finite(lost_bound_)) {
      res.bound = lost_bound_;
      return finish(Status::NumericalError);
    }
    res.bound = kInf;
    return finish(Status::Infeasible);
  }
  res.bound = std::min(bound, inc_obj_);
  const double gap_abs = inc_obj_ - res.bound;
  const bool closed = gap_abs <= params_.mip_abs_gap || gap_abs <= params_.mip_rel_gap * std::abs(inc_obj_) + 1e-9;
  return finish(closed ? Status::Optimal : Status::NumericalError);
}

}  // namespace

MiqpResult miqp_solve(const MipProblem& p, const SolverParams& params, const Logger& log, const Timer& timer,
                      const MiqpOptions& opt) {
  Miqp m(p, params, log, timer, opt);
  return m.run();
}

}  // namespace milp
