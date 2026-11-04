#include "milp/search/genres.hpp"

#include <algorithm>
#include <cmath>

#include "milp/cuts/cut.hpp"

namespace milp {

namespace {

constexpr double kBigBound = 1e9;  // larger bounds are treated as infinite (complementing with them is unsafe)
double conflict_tol(double rhs) { return 1e-6 * (1.0 + std::abs(rhs)); }  // as in conflict.cpp

// A linear constraint sum val x <= rhs (sparse, possibly with duplicate-free indices).
struct Lin {
  std::vector<Index> idx;
  std::vector<double> val;
  double rhs = 0.0;
  void clear() { idx.clear(), val.clear(), rhs = 0.0; }
  double coef(Index j) const {  // (a clause may have two literals on one column)
    double s = 0.0;
    for (size_t k = 0; k < idx.size(); ++k) {
      if (idx[k] == j) s += val[k];
    }
    return s;
  }
};

// Per-thread scratch indexed by column (sized once per problem size).
struct Scratch {
  std::vector<double> g;    // dense coefficients of the conflict constraint
  std::vector<char> in;     // column in the conflict's nonzero list
  std::vector<int> eff;     // per bound key: overriding trail entry (kUnset: trail.last)
  std::vector<Index> eff_touched;
  std::vector<double> d;    // dense accumulator for reason aggregation
  std::vector<char> din;
  void ensure(size_t n) {
    if (g.size() < n) {
      g.assign(n, 0.0), in.assign(n, 0), d.assign(n, 0.0), din.assign(n, 0);
      eff.assign(2 * n, kUnset);
    }
  }
  static constexpr int kUnset = -3;
};

class Resolver {
 public:
  Resolver(const BoundTrail& tr, const SparseMatrix& rows, const std::vector<double>& row_lb,
           const std::vector<double>& row_ub, const std::vector<char>& is_int, const std::vector<double>& glb,
           const std::vector<double>& gub, const GenresOptions& opt, GenresInfo& info)
      : tr_(tr), rows_(rows), row_lb_(row_lb), row_ub_(row_ub), is_int_(is_int), glb_(glb), gub_(gub), opt_(opt),
        info_(info) {
    thread_local Scratch s;
    s.ensure(is_int.size());
    s_ = &s;
  }
  ~Resolver() {
    for (Index j : nz_) s_->g[j] = 0.0, s_->in[j] = 0;
    for (Index k : s_->eff_touched) s_->eff[k] = Scratch::kUnset;
    s_->eff_touched.clear();
  }

  bool run(Conflict& out);

 private:
  // --- bounds ---
  double global(Index j, bool up) const {
    const double b = up ? gub_[j] : glb_[j];
    return std::abs(b) > kBigBound ? (b > 0 ? kInf : -kInf) : b;
  }
  int eff(Index j, bool up) const {
    const size_t key = 2 * static_cast<size_t>(j) + (up ? 1 : 0);
    const int k = s_->eff[key];
    return k == Scratch::kUnset ? tr_.last[key] : k;
  }
  void set_eff(Index j, bool up, int k) {
    const size_t key = 2 * static_cast<size_t>(j) + (up ? 1 : 0);
    if (s_->eff[key] == Scratch::kUnset) s_->eff_touched.push_back(static_cast<Index>(key));
    s_->eff[key] = k;
  }
  // Bound in effect: of the trail entry (no looser than the global bound).
  double bnd(Index j, bool up) const {
    const int k = eff(j, up);
    const double g = global(j, up);
    if (k < 0) return g;
    const double v = tr_.e[k].val;
    return up ? std::min(v, g) : std::max(v, g);
  }
  // Whether the bound in effect is tighter than the global one.
  bool local(Index j, bool up) const {
    const int k = eff(j, up);
    if (k < 0) return false;
    return up ? tr_.e[k].val < global(j, up) : tr_.e[k].val > global(j, up);
  }

  // --- the conflict constraint (dense g over nz_, rhs R_) ---
  void c_add(Index j, double v) {
    if (!s_->in[j]) s_->in[j] = 1, nz_.push_back(j);
    s_->g[j] += v;
  }
  void c_add(const Lin& l, double mult) {
    for (size_t k = 0; k < l.idx.size(); ++k) c_add(l.idx[k], mult * l.val[k]);
    R_ += mult * l.rhs;
    work_ += static_cast<long long>(l.idx.size());
  }
  // Drops exact zeros from nz_ and returns the number of nonzeros.
  size_t c_compact() {
    size_t w = 0;
    for (Index j : nz_) {
      if (s_->g[j] != 0.0) nz_[w++] = j;
      else s_->in[j] = 0;
    }
    nz_.resize(w);
    return w;
  }
  void c_save() {
    saved_.clear();
    for (Index j : nz_) {
      if (s_->g[j] != 0.0) saved_.idx.push_back(j), saved_.val.push_back(s_->g[j]);
    }
    saved_.rhs = R_;
  }
  void c_restore() {
    for (Index j : nz_) s_->g[j] = 0.0, s_->in[j] = 0;
    nz_.clear();
    R_ = 0.0;
    c_add(saved_, 1.0);
  }
  // Minimal activity of the conflict under the bounds in effect (false: unbounded).
  bool c_minact(double& act, double& absum) const {
    act = 0.0, absum = 0.0;
    for (Index j : nz_) {
      const double a = s_->g[j];
      if (a == 0.0) continue;
      const double b = bnd(j, a < 0);
      if (!is_finite(b)) return false;
      act += a * b;
      absum += std::abs(a * b);
    }
    return true;
  }
  // Infeasible under the bounds in effect, with a margin (so that the cleaned result still is).
  bool c_infeasible() const {
    double act, absum;
    if (!c_minact(act, absum)) return false;
    return act > R_ + conflict_tol(R_) + 1e-9 * absum;
  }
  void c_normalize() {
    double mx = 0.0;
    for (Index j : nz_) mx = std::max(mx, std::abs(s_->g[j]));
    if (!(mx > 0.0)) return;
    int e = 0;
    std::frexp(mx, &e);
    if (e == 0 || e == 1) return;
    for (Index j : nz_) s_->g[j] = std::ldexp(s_->g[j], -e);
    R_ = std::ldexp(R_, -e);
  }
  void coef_tightening();

  // --- reasons ---
  // The reason of trail entry k as sum val x <= rhs implying its bound; false if unavailable.
  bool reason(int k, Lin& r);
  bool conflict_as_lin(const Conflict& c, Lin& r);
  bool resolve(int k);
  bool reduce(const Lin& r, Index j, bool up, Lin& out);

  const BoundTrail& tr_;
  const SparseMatrix& rows_;
  const std::vector<double>&row_lb_, &row_ub_;
  const std::vector<char>& is_int_;
  const std::vector<double>&glb_, &gub_;
  const GenresOptions& opt_;
  GenresInfo& info_;
  Scratch* s_ = nullptr;
  std::vector<Index> nz_;
  double R_ = 0.0;
  double cutoff_ = kInf;
  long long work_ = 0;
  Lin saved_, reason_, reduced_, sub_;
};

bool Resolver::conflict_as_lin(const Conflict& c, Lin& r) {
  r.clear();
  cutoff_ = std::min(cutoff_, c.cutoff);
  if (!c.is_clause()) {
    r.idx = c.idx, r.val = c.val;
    r.rhs = c.rhs + conflict_tol(c.rhs);  // the slack the conflict propagates with
    return true;
  }
  // A clause "not all of x_k <= v_k (upper) / x_k >= v_k (lower)" whose literals all sit at global
  // bounds (v_k = glb_k for upper, gub_k for lower literals; binaries always) is linear: every term
  // of  sum_upper (x_k - glb_k) + sum_lower (gub_k - x_k) >= 1  is nonnegative, and the one of a
  // literal that does not hold is >= 1. As <=:  sum_lower x_k - sum_upper x_k <= sum_lower gub_k
  // - sum_upper glb_k - 1. Other literals (inside the global range) have no such form: fail.
  double rhs = -1.0;
  for (size_t k = 0; k < c.idx.size(); ++k) {
    const Index j = c.idx[k];
    if (!is_int_[j]) return false;
    const double gb = c.upper[k] ? glb_[j] : gub_[j];
    if (c.val[k] != gb || std::abs(gb) > kBigBound) return false;
    if (c.upper[k]) {
      r.idx.push_back(j), r.val.push_back(-1.0);
      rhs -= gb;
    } else {
      r.idx.push_back(j), r.val.push_back(1.0);
      rhs += gb;
    }
  }
  r.rhs = rhs;
  return true;
}

bool Resolver::reason(int k, Lin& r) {
  const BoundTrail::Entry& e = tr_.e[k];
  r.clear();
  if (e.kind == BoundTrail::kRow) {
    double a = 0.0;
    for (Index q = rows_.begin(e.row); q < rows_.end(e.row); ++q) {
      if (rows_.index[q] == e.col) a = rows_.value[q];
    }
    if (a == 0.0) return false;
    // x <= with a > 0, or x >= with a < 0: the row_ub side; else the row_lb side (negated).
    const bool upper_side = (a > 0) == (e.upper != 0);
    const double side = upper_side ? row_ub_[e.row] : -row_lb_[e.row];
    if (!is_finite(side) || std::abs(side) > 1e12) return false;
    const double sgn = upper_side ? 1.0 : -1.0;
    for (Index q = rows_.begin(e.row); q < rows_.end(e.row); ++q) {
      r.idx.push_back(rows_.index[q]), r.val.push_back(sgn * rows_.value[q]);
    }
    r.rhs = side;
    return true;
  }
  if (e.kind == BoundTrail::kConflict) {
    if (!conflict_as_lin(*e.conf, r)) {
      info_.reason = GenresInfo::kClauseReason;
      return false;
    }
    return true;
  }
  return false;
}

// Reduces reason r (implying the bound `up` of integer column j) so that its combination with the
// conflict stays infeasible: continuous columns whose bound in effect is local are resolved out
// with their reasons, the others are weakened to their global bounds, then a c-MIR rounding with
// divisor |a_j| is applied (integers complemented at the global bound their local bound sits at).
bool Resolver::reduce(const Lin& r, Index j, bool up, Lin& out) {
  // Dense aggregation of the reason.
  std::vector<double>& d = s_->d;
  std::vector<char>& din = s_->din;
  std::vector<Index> list;
  double rhs = r.rhs;
  auto add = [&](Index k, double v) {
    if (!din[k]) din[k] = 1, list.push_back(k);
    d[k] += v;
  };
  auto cleanup = [&]() {
    for (Index k : list) d[k] = 0.0, din[k] = 0;
  };
  for (size_t k = 0; k < r.idx.size(); ++k) add(r.idx[k], r.val[k]);
  // Resolve local continuous bounds out (latest first).
  for (int sub = 0;; ++sub) {
    int best = -1;
    Index bc = -1;
    for (Index k : list) {
      const double a = d[k];
      if (k == j || is_int_[k] || a == 0.0 || !local(k, a < 0)) continue;
      const int e = eff(k, a < 0);
      if (e > best) best = e, bc = k;
    }
    if (best < 0) break;
    if (sub >= 50 || work_ > opt_.max_work) {
      info_.reason = GenresInfo::kLimit;
      cleanup();
      return false;
    }
    if (tr_.e[best].kind == BoundTrail::kDecision || !reason(best, sub_)) {
      if (tr_.e[best].kind == BoundTrail::kDecision) info_.reason = GenresInfo::kContDecision;
      cleanup();
      return false;
    }
    const double ak = d[bc], bk = sub_.coef(bc);
    if (!(ak * bk < 0.0)) {
      cleanup();
      return false;
    }
    const double mult = std::abs(ak) / std::abs(bk);
    for (size_t k = 0; k < sub_.idx.size(); ++k) add(sub_.idx[k], mult * sub_.val[k]);
    rhs += mult * sub_.rhs;
    d[bc] = 0.0;
    work_ += static_cast<long long>(sub_.idx.size());
  }
  const double aj = d[j];
  if (aj == 0.0 || (up ? !(aj > 0) : !(aj < 0))) {
    cleanup();
    return false;
  }
  // Weaken continuous columns (their bounds in effect are global) and complement integers.
  const double delta = std::abs(aj);
  struct Term {
    Index k;
    double a;    // coefficient of z in the transformed constraint
    double b;    // complementing bound
    bool upper;  // x = b - z (else x = b + z)
  };
  std::vector<Term> terms;
  double beta = rhs;
  for (Index k : list) {
    const double a = d[k];
    if (a == 0.0) continue;
    if (!is_int_[k]) {
      const double gb = global(k, a < 0);
      if (!is_finite(gb)) {
        cleanup();
        return false;
      }
      beta -= a * gb;  // a x_k >= a gb
      continue;
    }
    // The bound the reason uses (lower for a > 0); complement at the global bound it equals,
    // otherwise (a general integer strictly inside its global range) at the nearer global bound.
    const bool use_up = a < 0;
    const double b = bnd(k, use_up), gl = global(k, false), gu = global(k, true);
    bool cu;
    if (k == j) {
      cu = !is_finite(gl);
    } else if (b == gl) {
      cu = false;
    } else if (b == gu) {
      cu = true;
    } else {
      cu = !is_finite(gl) || (is_finite(gu) && gu - b < b - gl);
    }
    const double cb = cu ? gu : gl;
    if (!is_finite(cb)) {
      cleanup();
      return false;
    }
    beta -= a * cb;
    terms.push_back({k, cu ? -a : a, cb, cu});
  }
  cleanup();
  // MIR on sum a'_t z_t <= beta (z >= 0 integer) with divisor delta:
  //   sum delta F(a'_t / delta) z_t <= delta floor(beta / delta),  F(u) = floor(u) + max(0, frac(u) - f0) / (1 - f0).
  // Without a usable fractional part (f0 ~ 0 or 1) the reason is kept as is (weakened only).
  out.clear();
  const double q = beta / delta;
  const double f0 = q - std::floor(q);
  const bool round = f0 > 1e-6 && f0 < 1.0 - 1e-6 && std::abs(q) < 1e9;
  if (round) ++info_.mir;
  double orhs = round ? delta * std::floor(q) : beta;
  for (const Term& t : terms) {
    double c = t.a;
    if (round) {
      const double u = t.a / delta;
      const double fu = u - std::floor(u);
      c = delta * (std::floor(u) + std::max(0.0, fu - f0) / (1.0 - f0));
    }
    if (c == 0.0) continue;
    // Back to x: z = x - b gives c x <= .. + c b; z = b - x gives -c x <= .. - c b.
    out.idx.push_back(t.k);
    out.val.push_back(t.upper ? -c : c);
    orhs += t.upper ? -c * t.b : c * t.b;
  }
  out.rhs = orhs;
  return std::isfinite(orhs);
}

void Resolver::coef_tightening() {
  double maxact = 0.0;
  for (Index k : nz_) {
    const double a = s_->g[k];
    if (a == 0.0) continue;
    const double b = global(k, a > 0);
    if (!is_finite(b)) return;
    maxact += a * b;
  }
  const double eps = 1e-9 * (1.0 + std::abs(R_));
  for (Index k : nz_) {
    const double a = s_->g[k];
    if (a == 0.0 || !is_int_[k]) continue;
    const double gl = global(k, false), gu = global(k, true);
    if (gu - gl != 1.0) continue;
    if (a > 0 && maxact - a < R_ - eps) {
      // x_k = gl makes the constraint redundant: lower a to the largest value that still matters.
      const double d = R_ - (maxact - a);
      s_->g[k] = a - d;
      R_ -= d * gu;
      maxact -= d * gu;
    } else if (a < 0 && maxact + a < R_ - eps) {
      const double d = R_ - (maxact + a);
      s_->g[k] = a + d;
      R_ += d * gl;
      maxact += d * gl;
    }
  }
}

bool Resolver::resolve(int k) {
  const BoundTrail::Entry& e = tr_.e[k];
  const Index j = e.col;
  const bool up = e.upper != 0;
  if (!reason(k, reason_)) return false;
  const double gj = s_->g[j], aj = reason_.coef(j);
  if (!(gj * aj < 0.0)) return false;
  c_save();
  set_eff(j, up, e.prev);
  // Plain combination (exact for continuous columns and for tight integer propagations).
  c_add(reason_, std::abs(gj) / std::abs(aj));
  s_->g[j] = 0.0;
  bool ok = c_infeasible();
  if (!ok && is_int_[j]) {
    c_restore();
    ++info_.reduced;
    if (reduce(reason_, j, up, reduced_)) {
      const double rj = reduced_.coef(j);
      if (rj * gj < 0.0) {
        c_add(reduced_, std::abs(gj) / std::abs(rj));
        s_->g[j] = 0.0;
        ok = c_infeasible();
      }
    }
    if (!ok && info_.reason == GenresInfo::kNone) info_.reason = GenresInfo::kReduction;
  } else if (!ok) {
    info_.reason = GenresInfo::kReduction;
  }
  if (!ok) {
    c_restore();
    set_eff(j, up, k);
    return false;
  }
  if (static_cast<long long>(c_compact()) > 2LL * std::max<Index>(opt_.max_nnz, 1)) {
    info_.reason = GenresInfo::kLimit;
    c_restore();
    set_eff(j, up, k);
    return false;
  }
  c_normalize();
  if (opt_.coef_tightening) {
    c_save();
    coef_tightening();
    if (!c_infeasible()) c_restore();
  }
  return true;
}

bool Resolver::run(Conflict& out) {
  // The failing constraint.
  switch (tr_.fail) {
    case BoundTrail::kFailRow: {
      const Index i = tr_.fail_row;
      const double side = tr_.fail_upper ? row_ub_[i] : -row_lb_[i];
      if (!is_finite(side) || std::abs(side) > 1e12) break;
      const double sgn = tr_.fail_upper ? 1.0 : -1.0;
      for (Index q = rows_.begin(i); q < rows_.end(i); ++q) c_add(rows_.index[q], sgn * rows_.value[q]);
      R_ = side;
      break;
    }
    case BoundTrail::kFailConflict:
      if (conflict_as_lin(*tr_.fail_conf, reason_)) c_add(reason_, 1.0);
      break;
    case BoundTrail::kFailCrossing: {
      // The later of the two crossing bounds (not a decision): its reason is infeasible without it.
      const Index j = tr_.fail_col;
      const int kl = tr_.last[2 * static_cast<size_t>(j)], ku = tr_.last[2 * static_cast<size_t>(j) + 1];
      int k = std::max(kl, ku);
      if (k >= 0 && tr_.e[k].kind == BoundTrail::kDecision) k = std::min(kl, ku);
      if (k < 0 || tr_.e[k].kind == BoundTrail::kDecision || !reason(k, reason_)) break;
      c_add(reason_, 1.0);
      set_eff(j, tr_.e[k].upper != 0, tr_.e[k].prev);
      break;
    }
    default:
      break;
  }
  if (nz_.empty()) {
    if (info_.reason == GenresInfo::kNone) info_.reason = GenresInfo::kNoStart;
    return false;
  }
  c_compact();
  c_normalize();
  if (!c_infeasible()) {
    info_.reason = GenresInfo::kNoStart;
    return false;
  }
  // Cheap exit for the common case (LP proofs on many models) that every local bound the constraint
  // uses is a decision: nothing to resolve. (Decisions precede the propagated bounds on the trail, so
  // relaxing a decision exposes only earlier decisions or the global bound.)
  bool resolvable = false;
  for (Index j : nz_) {
    const double a = s_->g[j];
    if (a != 0.0 && local(j, a < 0) && tr_.e[eff(j, a < 0)].kind != BoundTrail::kDecision) {
      resolvable = true;
      break;
    }
  }
  if (!resolvable) {
    info_.reason = GenresInfo::kNoStep;
    return false;
  }
  GenresInfo::Outcome outcome = GenresInfo::kFail;
  for (;;) {
    // The latest local bound the conflict uses.
    int best = -1, maxlv = -1, atmax = 0;
    for (Index j : nz_) {
      const double a = s_->g[j];
      if (a == 0.0 || !local(j, a < 0)) continue;
      const int k = eff(j, a < 0);
      best = std::max(best, k);
      const int lv = tr_.e[k].level;
      if (lv > maxlv) maxlv = lv, atmax = 1;
      else if (lv == maxlv) ++atmax;
    }
    work_ += static_cast<long long>(nz_.size());
    if (work_ > 4 * opt_.max_work) {
      info_.reason = GenresInfo::kLimit;
      if (info_.steps == 0 || !opt_.partial) return false;
      outcome = GenresInfo::kPartial;
      break;
    }
    if (best < 0) {
      outcome = GenresInfo::kGlobal;  // infeasible under the global bounds
      break;
    }
    const BoundTrail::Entry& e = tr_.e[best];
    // Not needed: relax it to the previous bound.
    set_eff(e.col, e.upper != 0, e.prev);
    if (c_infeasible()) continue;
    set_eff(e.col, e.upper != 0, best);
    if (info_.steps > 0 && atmax == 1) {
      outcome = GenresInfo::kUip;
      break;
    }
    if (e.kind == BoundTrail::kDecision) {
      if (info_.steps == 0) {
        info_.reason = GenresInfo::kNoStep;
        return false;
      }
      outcome = GenresInfo::kDecision;
      break;
    }
    if (info_.steps >= opt_.max_steps || work_ > opt_.max_work) {
      info_.reason = GenresInfo::kLimit;
      if (info_.steps == 0 || !opt_.partial) return false;
      outcome = GenresInfo::kPartial;
      break;
    }
    if (!resolve(best)) {
      if (info_.reason == GenresInfo::kNone) info_.reason = GenresInfo::kNumerics;
      if (info_.steps == 0 || !opt_.partial) return false;
      outcome = GenresInfo::kPartial;
      break;
    }
    ++info_.steps;
  }
  if (info_.steps == 0) {
    info_.reason = GenresInfo::kNoStep;
    return false;
  }
  // Clean against the global bounds and check that it proves the node infeasible.
  Cut cc;
  for (Index j : nz_) {
    if (s_->g[j] != 0.0) cc.idx.push_back(j), cc.val.push_back(s_->g[j]);
  }
  cc.rhs = R_;
  if (!clean_cut(cc, glb_, gub_, opt_.max_nnz)) {
    info_.reason = GenresInfo::kLimit;
    return false;
  }
  double act = 0.0;
  for (size_t k = 0; k < cc.idx.size(); ++k) {
    const Index j = cc.idx[k];
    const bool up = cc.val[k] < 0;
    const int t = tr_.last[2 * static_cast<size_t>(j) + (up ? 1 : 0)];
    double b = global(j, up);
    if (t >= 0) b = up ? std::min(b, tr_.e[t].val) : std::max(b, tr_.e[t].val);
    if (!is_finite(b)) {
      info_.reason = GenresInfo::kFinalCheck;
      return false;
    }
    act += cc.val[k] * b;
  }
  if (!(act > cc.rhs + conflict_tol(cc.rhs))) {
    info_.reason = GenresInfo::kFinalCheck;
    return false;
  }
  out.idx = std::move(cc.idx);
  out.val = std::move(cc.val);
  out.upper.clear();
  out.rhs = cc.rhs;
  out.cutoff = cutoff_;
  info_.outcome = outcome;
  return true;
}

}  // namespace

bool analyze_genres(const BoundTrail& trail, const SparseMatrix& rows, const std::vector<double>& row_lb,
                    const std::vector<double>& row_ub, const std::vector<char>& is_int, const std::vector<double>& glb,
                    const std::vector<double>& gub, const GenresOptions& opt, Conflict& out, GenresInfo* info) {
  GenresInfo local_info;
  GenresInfo& inf = info ? *info : local_info;
  inf = GenresInfo();
  if (!trail.usable()) {
    inf.reason = GenresInfo::kNoStart;
    return false;
  }
  Resolver r(trail, rows, row_lb, row_ub, is_int, glb, gub, opt, inf);
  return r.run(out);
}

}  // namespace milp
