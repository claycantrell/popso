// Tests for conflict analysis: Farkas proofs from the dual simplex, dual proof
// constraints (validity against every feasible mixed-integer point better than
// the cutoff), the conflict pool's propagation, and end-to-end solves.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <random>
#include <string>

#include "milp/core/model.hpp"
#include "milp/lp/simplex.hpp"
#include "milp/search/conflict.hpp"
#include "milp/search/genres.hpp"
#include "milp/search/trail.hpp"
#include "milp/search/mip_problem.hpp"
#include "milp/search/propagation.hpp"
#include "milp/solver/solver.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

int scaled(int base) {
#ifdef _MSC_VER
#pragma warning(suppress : 4996)
#endif
  const char* v = std::getenv("MILP_TEST_SCALE");
  return v ? base * std::max(1, std::atoi(v)) : base;
}

void set_env(const char* name, bool on) {
#ifdef _WIN32
  _putenv_s(name, on ? "1" : "");
#else
  if (on) setenv(name, "1", 1);
  else unsetenv(name);
#endif
}

// Random MIP: integers in [-1, 3] or binaries, bounded continuous variables,
// random rows around a random point (mostly feasible).
Model random_mip(std::mt19937& rng, int ni, int nc, int nr) {
  ModelBuilder b;
  std::uniform_int_distribution<int> coef(-9, 9);
  std::vector<double> x0;
  for (int j = 0; j < ni; ++j) {
    const bool bin = rng() % 2;
    const double lo = bin ? 0 : -1.0 + static_cast<int>(rng() % 2), up = bin ? 1 : 3;
    b.add_col("i" + std::to_string(j), coef(rng), lo, up, bin ? VarType::Binary : VarType::Integer);
    x0.push_back(lo + static_cast<int>(rng() % static_cast<unsigned>(up - lo + 1)));
  }
  for (int j = 0; j < nc; ++j) {
    const double lo = rng() % 3 == 0 ? -2.0 : 0.0, up = 4.0 + rng() % 6;
    b.add_col("c" + std::to_string(j), coef(rng) * 0.5, lo, up);
    x0.push_back(lo + (up - lo) * (rng() % 100) / 100.0);
  }
  for (int r = 0; r < nr; ++r) {
    std::vector<std::pair<Index, double>> terms;
    for (int j = 0; j < ni + nc; ++j) {
      if (rng() % 3 == 0) continue;
      double a = coef(rng);
      if (a == 0.0) a = 1.0;
      if (j >= ni) a += (rng() % 10) / 10.0;
      terms.emplace_back(j, a);
    }
    if (terms.empty()) continue;
    double act = 0.0;
    for (auto& t : terms) act += t.second * x0[t.first];
    const double slack = (rng() % 5) * 0.5;
    const int kind = rng() % 3;
    const double lo = kind == 0 ? -kInf : act - slack - (kind == 2 ? rng() % 3 : 0);
    const Index row = b.add_row("r" + std::to_string(r), lo, kind == 1 ? kInf : act + slack);
    for (auto& t : terms) b.add_coef(row, t.first, t.second);
  }
  return b.build();
}

// Max of sum val x over points with the integers fixed at x_int, satisfying the
// rows and, if finite, cost x <= cutoff. Returns -inf if there is no such point.
double max_lhs(const MipProblem& p, const std::vector<double>& x_int, const Conflict& c, double cutoff) {
  std::vector<double> lb = p.col_lb, ub = p.col_ub, obj(p.num_cols(), 0.0);
  for (Index j = 0; j < p.num_cols(); ++j) {
    if (p.is_int[j]) lb[j] = ub[j] = x_int[j];
  }
  for (size_t k = 0; k < c.idx.size(); ++k) obj[c.idx[k]] -= c.val[k];
  SparseMatrix A = p.A;
  std::vector<double> rlb = p.row_lb, rub = p.row_ub;
  if (is_finite(cutoff)) {
    // Append the objective row as an extra row: rebuild via triplets.
    std::vector<Index> maj, mnr;
    std::vector<double> val;
    for (Index j = 0; j < p.num_cols(); ++j) {
      for (Index k = p.A.begin(j); k < p.A.end(j); ++k) maj.push_back(j), mnr.push_back(p.A.index[k]), val.push_back(p.A.value[k]);
      if (p.cost[j] != 0.0) maj.push_back(j), mnr.push_back(p.num_rows()), val.push_back(p.cost[j]);
    }
    A = SparseMatrix::from_triplets(p.num_cols(), p.num_rows() + 1, maj, mnr, val);
    rlb.push_back(-kInf);
    rub.push_back(cutoff);
  }
  DualSimplex lp;
  lp.load(A, obj, lb, ub, rlb, rub);
  LpOptions o;
  o.perturb = false;
  const LpStatus st = lp.solve(o);
  if (st == LpStatus::Infeasible) return -kInf;
  if (st != LpStatus::Optimal) return kInf;
  return -lp.objective();
}

template <class F>
void enumerate_ints(const MipProblem& p, F&& f) {
  const Index n = p.num_cols();
  std::vector<double> x(n, 0.0);
  for (Index j = 0; j < n; ++j) x[j] = p.is_int[j] ? p.col_lb[j] : 0.0;
  for (;;) {
    f(x);
    Index j = 0;
    for (; j < n; ++j) {
      if (!p.is_int[j]) continue;
      if (x[j] < p.col_ub[j]) {
        x[j] += 1.0;
        break;
      }
      x[j] = p.col_lb[j];
    }
    if (j == n) return;
  }
}

}  // namespace

TEST(conflict_farkas_proof_simple) {
  // x + y >= 3 with x, y in [0, 1] is infeasible.
  SparseMatrix A = SparseMatrix::from_triplets(2, 1, {0, 1}, {0, 0}, {1.0, 1.0});
  DualSimplex lp;
  lp.load(A, {1.0, 1.0}, {0.0, 0.0}, {1.0, 1.0}, {3.0}, {kInf});
  LpOptions o;
  REQUIRE(lp.solve(o) == LpStatus::Infeasible);
  std::vector<double> y;
  REQUIRE(lp.farkas_proof(y));
  REQUIRE(y.size() == 1u);
  CHECK(std::abs(y[0]) > 0.0);
  Conflict c;
  const std::vector<char> is_int{1, 1};
  REQUIRE(infeasibility_conflict(lp, y, {1.0, 1.0}, is_int, {0.0, 0.0}, {1.0, 1.0}, ConflictOptions{}, c));
  // The conflict is -x - y <= -3 up to scaling: infeasible over [0, 1]^2.
  std::vector<double> lb{0.0, 0.0}, ub{1.0, 1.0};
  CHECK(!propagate_conflict(c, is_int, lb, ub, nullptr, nullptr));
  // An optimal solve has no proof.
  lp.set_row_bounds(0, 1.0, kInf);
  REQUIRE(lp.solve(o) == LpStatus::Optimal);
  CHECK(!lp.farkas_proof(y));
}

TEST(conflict_dual_proofs_valid_random) {
  std::mt19937 rng(5);
  int infeas = 0, cutoff = 0, bad = 0;
  for (int t = 0; t < scaled(400); ++t) {
    Model m = random_mip(rng, 3 + rng() % 4, rng() % 3, 2 + rng() % 4);
    MipProblem p;
    if (!MipProblem::from_model(m, &p)) continue;
    DualSimplex lp;
    lp.load(p.A, p.cost, p.col_lb, p.col_ub, p.row_lb, p.row_ub);
    LpOptions o;
    o.perturb = t % 2 == 0;
    const LpStatus st0 = lp.solve(o);
    if (st0 != LpStatus::Optimal) continue;
    const double root = lp.objective();
    // Random local bounds (branching on integers) and a cutoff a bit above the root bound.
    for (int k = 0; k < 6; ++k) {
      const Index j = static_cast<Index>(rng() % p.num_cols());
      if (!p.is_int[j]) continue;
      const double l = lp.col_lower(j), u = lp.col_upper(j);
      const double v = l + static_cast<int>(rng() % static_cast<unsigned>(u - l + 1));
      if (rng() % 2) lp.set_col_bounds(j, v, u);
      else lp.set_col_bounds(j, l, v);
    }
    const double co = root + (rng() % 20) * 0.5;
    o.objective_cutoff = rng() % 2 ? co : kInf;
    const LpStatus st = lp.solve(o);
    Conflict c;
    bool ok = false;
    double used_cutoff = kInf;
    if (st == LpStatus::Infeasible) {
      std::vector<double> y;
      REQUIRE(lp.farkas_proof(y));
      ok = infeasibility_conflict(lp, y, p.cost, p.is_int, p.col_lb, p.col_ub, ConflictOptions{}, c);
      infeas += ok;
    } else if (st == LpStatus::CutoffReached || (st == LpStatus::Optimal && lp.objective() >= co)) {
      std::vector<double> mu = lp.row_duals();
      for (double& v : mu) v = -v;
      ConflictOptions copt;
      copt.relax_ints = rng() % 2;
      ok = derive_conflict(lp, mu, 1.0, p.cost, co, p.is_int, p.col_lb, p.col_ub, copt, c);
      used_cutoff = co;
      cutoff += ok;
    }
    if (!ok) continue;
    // The conflict prunes the node ...
    std::vector<double> lb(p.num_cols()), ub(p.num_cols());
    for (Index j = 0; j < p.num_cols(); ++j) lb[j] = lp.col_lower(j), ub[j] = lp.col_upper(j);
    CHECK(!propagate_conflict(c, p.is_int, lb, ub, nullptr, nullptr));
    // ... and holds at every feasible point (better than the cutoff).
    double worst = -kInf;
    enumerate_ints(p, [&](const std::vector<double>& x) { worst = std::max(worst, max_lhs(p, x, c, used_cutoff)); });
    if (worst > c.rhs + 1e-6 * (1.0 + std::abs(c.rhs))) {
      ++bad;
      std::fprintf(stderr, "  invalid conflict (%s): max lhs %g > rhs %g\n", is_finite(used_cutoff) ? "cutoff" : "infeasible",
                   worst, c.rhs);
    }
  }
  std::printf("  %d infeasibility and %d cutoff proofs checked\n", infeas, cutoff);
  CHECK(infeas > 20);
  CHECK(cutoff > 20);
  CHECK_EQ(bad, 0);
}

TEST(conflict_pool_propagates_and_ages) {
  const std::vector<char> is_int(4, 1);
  const std::vector<double> glb(4, 0.0), gub(4, 1.0);
  ConflictPool pool(4, 2);
  auto c1 = std::make_shared<Conflict>();
  c1->idx = {0, 1, 2};  // x0 + x1 + x2 <= 1
  c1->val = {1.0, 1.0, 1.0};
  c1->rhs = 1.0;
  pool.add(c1, 0);
  auto c2 = std::make_shared<Conflict>();
  c2->idx = {2, 3};  // x2 - x3 <= 0
  c2->val = {1.0, -1.0};
  c2->rhs = 0.0;
  pool.add(c2, 0);
  std::vector<double> lb = glb, ub = gub;
  std::vector<Index> t;
  lb[0] = 1.0;
  REQUIRE(pool.propagate(is_int, glb, gub, lb, ub, {0}, t, 1));
  CHECK_EQ(ub[1], 0.0);
  CHECK_EQ(ub[2], 0.0);
  CHECK_EQ(lb[3], 0.0);
  // x2 = 1 and x0 = 1 contradict c1.
  lb = glb, ub = gub;
  lb[0] = lb[2] = 1.0;
  t.clear();
  CHECK(!pool.propagate(is_int, glb, gub, lb, ub, {0, 2}, t, 2));
  // Propagation through the chain: x2 = 1 forces x3 = 1.
  lb = glb, ub = gub;
  lb[2] = 1.0;
  t.clear();
  REQUIRE(pool.propagate(is_int, glb, gub, lb, ub, {2}, t, 3));
  CHECK_EQ(lb[3], 1.0);
  CHECK_EQ(ub[0], 0.0);
  // Over capacity: the least recently useful conflict goes first.
  auto c3 = std::make_shared<Conflict>();
  c3->idx = {3};
  c3->val = {1.0};
  c3->rhs = 5.0;  // never propagates
  pool.add(c3, 4);
  CHECK(pool.size() <= 2u);
  lb = glb, ub = gub;
  lb[2] = 1.0;
  t.clear();
  REQUIRE(pool.propagate(is_int, glb, gub, lb, ub, {2}, t, 5));
  CHECK_EQ(lb[3], 1.0);  // c2 (used at node 3) survived
}

TEST(conflict_pool_keeps_trail_reasons_alive) {
  // A conflict that implied a bound on the trail stays alive while the trail refers to it, also when adding
  // conflicts shrinks the pool and drops it (neos-2657525-crna: the propagation analysis read a freed conflict).
  const std::vector<char> is_int(4, 1);
  const std::vector<double> glb(4, 0.0), gub(4, 1.0);
  ConflictPool pool(4, 2);
  std::weak_ptr<const Conflict> weak;
  {
    auto c = std::make_shared<Conflict>();
    c->idx = {0, 1};  // x0 + x1 <= 1
    c->val = {1.0, 1.0};
    c->rhs = 1.0;
    weak = c;
    pool.add(c, 0);
  }
  BoundTrail trail;
  trail.init(4, 100);
  std::vector<double> lb = glb, ub = gub;
  std::vector<Index> t;
  lb[0] = 1.0;
  trail.push(0, false, 1.0, BoundTrail::kDecision);
  REQUIRE(pool.propagate(is_int, glb, gub, lb, ub, {0}, t, 1, &trail));
  CHECK_EQ(ub[1], 0.0);
  REQUIRE(trail.e.size() == 2 && trail.e[1].kind == BoundTrail::kConflict);
  for (int k = 0; k < 4; ++k) {  // newer conflicts: the first one is now the least recently used
    auto d = std::make_shared<Conflict>();
    d->idx = {2, 3};
    d->val = {1.0, 1.0};
    d->rhs = 1.0;
    pool.add(d, 10 + k);
  }
  CHECK(!weak.expired());
  CHECK(trail.e[1].conf->idx.size() == 2);
  trail.clear();
  CHECK(weak.expired());  // the pool had dropped it
}

TEST(conflict_end_to_end_same_optimum) {
  // Solves with and without conflict analysis must agree (opportunistic and deterministic).
  std::mt19937 rng(77);
  SolverParams params;
  params.verbosity = 0;
  params.mip_rel_gap = 1e-9;
  params.mip_abs_gap = 1e-9;
  params.node_limit = 20000;
  params.presolve = false;
  int compared = 0;
  for (int t = 0; t < scaled(40); ++t) {
    Model m = random_mip(rng, 25 + rng() % 15, rng() % 6, 10 + rng() % 8);
    params.threads = 1 + t % 2;
    params.deterministic = t % 4 == 1;
    set_env("MILP_NO_CONFLICT", true);
    SolveResult r0 = solve(m, params);
    set_env("MILP_NO_CONFLICT", false);
    SolveResult r1 = solve(m, params);
    if (r0.status == Status::Infeasible) {
      CHECK(r1.status == Status::Infeasible);
      continue;
    }
    if (r0.status != Status::Optimal || r1.status != Status::Optimal) continue;
    CHECK_NEAR(r1.objective, r0.objective, 1e-6 * std::max(1.0, std::abs(r0.objective)));
    ++compared;
  }
  std::printf("  %d instances compared\n", compared);
  CHECK(compared > 10);
}

TEST(conflict_propagation_clauses_valid_by_enumeration) {
  // Random pure-integer problems; random dives (one decision per level, row
  // and conflict-pool propagation after each). When propagation fails, the
  // UIP clauses from the trail must hold for every integer point satisfying
  // the rows, checked by enumeration. Learned clauses (and the rows as linear
  // conflicts) go into the pool, so later clauses also have conflict reasons.
  std::mt19937 rng(11);
  int clauses = 0, from_pool = 0, bad = 0, short_ones = 0;
  for (int t = 0; t < scaled(300); ++t) {
    Model m = random_mip(rng, 4 + rng() % 3, 0, 3 + rng() % 4);
    MipProblem p;
    if (!MipProblem::from_model(m, &p)) continue;
    const Index n = p.num_cols();
    std::vector<std::vector<double>> feasible;
    enumerate_ints(p, [&](const std::vector<double>& x) {
      std::vector<double> act(p.num_rows(), 0.0);
      for (Index j = 0; j < n; ++j) {
        for (Index k = p.A.begin(j); k < p.A.end(j); ++k) act[p.A.index[k]] += p.A.value[k] * x[j];
      }
      for (Index i = 0; i < p.num_rows(); ++i) {
        if (act[i] < p.row_lb[i] - 1e-9 || act[i] > p.row_ub[i] + 1e-9) return;
      }
      feasible.push_back(x);
    });
    DomainPropagator prop(p);
    BoundTrail trail;
    trail.init(n, 10000);
    ConflictPool pool(n, 200);
    const SparseMatrix rows = p.A.transpose();
    for (Index i = 0; i < p.num_rows(); ++i) {
      if (!is_finite(p.row_ub[i]) || rng() % 2) continue;
      auto c = std::make_shared<Conflict>();
      for (Index k = rows.begin(i); k < rows.end(i); ++k) c->idx.push_back(rows.index[k]), c->val.push_back(rows.value[k]);
      c->rhs = p.row_ub[i];
      pool.add(c, 0);
    }
    for (int dive = 0; dive < 25; ++dive) {
      std::vector<double> lb = p.col_lb, ub = p.col_ub;
      trail.clear();
      prop.trail = &trail;
      bool ok = true;
      // Rows and pool alternately from the seeds (a few rounds).
      auto propagate = [&](std::vector<Index> seeds) {
        std::vector<Index> tight, ctight;
        for (int pass = 0; pass < 4 && ok && !seeds.empty(); ++pass) {
          tight.clear();
          ok = prop.propagate(lb, ub, seeds, tight);
          if (!ok) break;
          seeds.insert(seeds.end(), tight.begin(), tight.end());
          ctight.clear();
          ok = pool.propagate(p.is_int, p.col_lb, p.col_ub, lb, ub, seeds, ctight, dive, &trail);
          seeds = ctight;
        }
      };
      // Even dives propagate after each decision; odd ones apply all decisions
      // first (as the search does with a node's path), so implied bounds get
      // the highest level among their reasons.
      const bool bulk = dive % 2 == 1;
      std::vector<Index> decided;
      for (int level = 0; level < 8 && ok; ++level) {
        std::vector<Index> open;
        for (Index j = 0; j < n; ++j) {
          if (lb[j] < ub[j]) open.push_back(j);
        }
        if (open.empty()) break;
        const Index j = open[rng() % open.size()];
        const double v = lb[j] + static_cast<int>(rng() % static_cast<unsigned>(ub[j] - lb[j] + 1));
        ++trail.level;
        if (rng() % 2) {
          if (v == ub[j]) continue;
          ub[j] = v;
          trail.push(j, true, v, BoundTrail::kDecision);
        } else {
          if (v == lb[j]) continue;
          lb[j] = v;
          trail.push(j, false, v, BoundTrail::kDecision);
        }
        if (bulk) decided.push_back(j);
        else propagate({j});
      }
      if (bulk && ok) propagate(decided);
      prop.trail = nullptr;
      if (ok) continue;
      if (trail.fail == BoundTrail::kFailConflict) ++from_pool;
      std::vector<Conflict> out;
      analyze_propagation_conflict(trail, rows, p.row_lb, p.row_ub, p.is_int, p.col_lb, p.col_ub, 50, 3, out);
      for (Conflict& c : out) {
        ++clauses;
        if (c.idx.size() <= 2) ++short_ones;
        for (const auto& x : feasible) {
          if (c.violated_by(x, 0.0)) {
            ++bad;
            break;
          }
        }
        pool.add(std::make_shared<Conflict>(std::move(c)), dive);
      }
    }
  }
  std::printf("  %d clauses checked (%d short, %d failures in the pool), %d invalid\n", clauses, short_ones, from_pool,
              bad);
  CHECK(clauses > 50);
  CHECK(from_pool > 5);
  CHECK_EQ(bad, 0);
}

TEST(conflict_clause_propagation) {
  // Clause: not (x0 >= 2 and x1 <= 0 and x2 >= 1).
  const std::vector<char> is_int(3, 1);
  const std::vector<double> glb(3, 0.0), gub(3, 3.0);
  auto c = std::make_shared<Conflict>();
  c->idx = {0, 1, 2};
  c->val = {2.0, 0.0, 1.0};
  c->upper = {0, 1, 0};
  ConflictPool pool(3, 10);
  pool.add(c, 0);
  std::vector<double> lb = glb, ub = gub;
  std::vector<Index> t;
  lb[0] = 2.0;
  REQUIRE(pool.propagate(is_int, glb, gub, lb, ub, {0}, t, 1));
  CHECK_EQ(ub[2], 3.0);  // two literals open: nothing
  ub[1] = 0.0;
  t.clear();
  REQUIRE(pool.propagate(is_int, glb, gub, lb, ub, {1}, t, 2));
  CHECK_EQ(ub[2], 0.0);  // unit: x2 <= 0
  // At another node (watches moved): all literals hold is infeasible.
  lb = glb, ub = gub;
  lb[2] = 1.0, ub[1] = 0.0;
  t.clear();
  REQUIRE(pool.propagate(is_int, glb, gub, lb, ub, {2, 1}, t, 3));
  CHECK_EQ(ub[0], 1.0);
  lb[0] = 2.0, ub[0] = 3.0;
  t.clear();
  CHECK(!pool.propagate(is_int, glb, gub, lb, ub, {0}, t, 4));
  std::vector<double> l2 = glb, u2 = gub;
  l2[0] = 3.0, u2[1] = 0.0;
  CHECK(propagate_conflict(*c, is_int, l2, u2, nullptr, nullptr));
  CHECK_EQ(u2[2], 0.0);
}

TEST(node_lp_failures_keep_subtrees) {
  // Node LPs that fail and cannot be recovered (forced with MILP_DEBUG_LP_FAIL: every 3rd node) are
  // branched on without their LP instead of dropped: same optimum, same status.
  std::mt19937 rng(91);
  SolverParams params;
  params.verbosity = 0;
  params.mip_rel_gap = 1e-9;
  params.mip_abs_gap = 1e-9;
  params.node_limit = 50000;
  params.presolve = false;
  int compared = 0, branched = 0;
  for (int t = 0; t < scaled(30); ++t) {
    Model m = random_mip(rng, 25 + rng() % 15, rng() % 6, 10 + rng() % 8);
    params.threads = 1 + t % 2;
    params.deterministic = t % 4 == 1;
    SolveResult r0 = solve(m, params);
#ifdef _WIN32
    _putenv_s("MILP_DEBUG_LP_FAIL", "3");
#else
    setenv("MILP_DEBUG_LP_FAIL", "3", 1);
#endif
    SolveResult r1 = solve(m, params);
#ifdef _WIN32
    _putenv_s("MILP_DEBUG_LP_FAIL", "");
#else
    unsetenv("MILP_DEBUG_LP_FAIL");
#endif
    // (large trees grow further without a third of their LPs: may run into the node limit)
    if ((r0.status != Status::Optimal && r0.status != Status::Infeasible) || r0.stats.nodes > 5000) continue;
    if (r1.status != r0.status) std::printf("  t %d: status %d vs %d, obj %g vs %g, nodes %lld vs %lld\n", t, (int)r0.status, (int)r1.status, r0.objective, r1.objective, r0.stats.nodes, r1.stats.nodes);
    CHECK(r1.status == r0.status);
    if (r0.status == Status::Optimal && r1.status == Status::Optimal) {
      CHECK_NEAR(r1.objective, r0.objective, 1e-6 * std::max(1.0, std::abs(r0.objective)));
      ++compared;
      branched += r1.stats.nodes > 3;
    }
  }
  std::printf("  %d instances compared (%d with more than 3 nodes)\n", compared, branched);
  CHECK(compared > 10);
  CHECK(branched > 3);
}

TEST(points_at_lp_tolerance_are_polished) {
  // LP points feasible only within the LP's primal tolerance miss the 1e-6 model check (timtab1: basic
  // row activities 1e-6 + rounding beyond their sides; the node was then branched on without an LP down
  // to all integers fixed and dropped, NUMERIC with a worse incumbent). MILP_DEBUG_POINT_SLACK moves the
  // first continuous column at a bound of every point to verify 1.01e-6 outside it: the continuous part is
  // re-solved with the integers fixed, so the optimum and status stay the same, and the answer is
  // feasible.
  std::mt19937 rng(17);
  SolverParams params;
  params.verbosity = 0;
  params.mip_rel_gap = 1e-9;
  params.mip_abs_gap = 1e-9;
  params.node_limit = 50000;
  int compared = 0;
  for (int t = 0; t < scaled(24); ++t) {
    Model m = random_mip(rng, 15 + rng() % 10, 3 + rng() % 4, 8 + rng() % 6);
    params.threads = 1 + t % 2;
    params.deterministic = t % 4 == 1;
    params.presolve = t % 3 != 0;
    SolveResult r0 = solve(m, params);
#ifdef _WIN32
    _putenv_s("MILP_DEBUG_POINT_SLACK", "1.01e-6");
#else
    setenv("MILP_DEBUG_POINT_SLACK", "1.01e-6", 1);
#endif
    SolveResult r1 = solve(m, params);
#ifdef _WIN32
    _putenv_s("MILP_DEBUG_POINT_SLACK", "");
#else
    unsetenv("MILP_DEBUG_POINT_SLACK");
#endif
    if (r0.status != Status::Optimal && r0.status != Status::Infeasible) continue;
    if (r1.status != r0.status) std::printf("  t %d: status %d vs %d, obj %g vs %g\n", t, (int)r0.status, (int)r1.status, r0.objective, r1.objective);
    CHECK(r1.status == r0.status);
    if (r0.status == Status::Optimal && r1.status == Status::Optimal) {
      CHECK_NEAR(r1.objective, r0.objective, 1e-6 * std::max(1.0, std::abs(r0.objective)));
      MipProblem p;
      REQUIRE(MipProblem::from_model(m, &p));
      CHECK(p.max_violation(r1.x) <= 1e-6);
      ++compared;
    }
  }
  std::printf("  %d instances compared\n", compared);
  CHECK(compared > 8);
}

namespace {

// Random dives as in conflict_propagation_clauses_valid_by_enumeration, analysed by cut-based
// conflict analysis. Pure integer problems (nc = 0) are checked by enumeration, mixed ones by an LP
// over the continuous columns per integer point; with an objective cutoff (at the optimum) only
// the points with cost <= cutoff count, so every optimal solution must satisfy every learned
// constraint. Learned constraints and clauses go into the pool and become reasons later.
struct GenresStats {
  int learned = 0, uip = 0, partial = 0, steps = 0, reduced = 0, mir = 0, bad = 0, calls = 0, nnz = 0;
  int fails[9] = {};
};

void genres_dives(std::mt19937& rng, int trials, int nc_max, bool cutoff_on, GenresStats& st) {
  for (int t = 0; t < trials; ++t) {
    const int nc = nc_max > 0 ? 1 + static_cast<int>(rng() % static_cast<unsigned>(nc_max)) : 0;
    Model m = random_mip(rng, 4 + rng() % 3, nc, 3 + rng() % 4);
    MipProblem p;
    if (!MipProblem::from_model(m, &p)) continue;
    const Index n = p.num_cols();
    // Integer points with a feasible continuous part, and the optimum.
    std::vector<std::vector<double>> pts;
    double opt = kInf;
    Conflict negc;  // max of -cost x over the continuous part
    for (Index j = 0; j < n; ++j) {
      if (p.cost[j] != 0.0) negc.idx.push_back(j), negc.val.push_back(-p.cost[j]);
    }
    enumerate_ints(p, [&](const std::vector<double>& x) {
      if (nc == 0) {
        std::vector<double> act(p.num_rows(), 0.0);
        for (Index j = 0; j < n; ++j) {
          for (Index k = p.A.begin(j); k < p.A.end(j); ++k) act[p.A.index[k]] += p.A.value[k] * x[j];
        }
        for (Index i = 0; i < p.num_rows(); ++i) {
          if (act[i] < p.row_lb[i] - 1e-9 || act[i] > p.row_ub[i] + 1e-9) return;
        }
        double obj = 0.0;
        for (Index j = 0; j < n; ++j) obj += p.cost[j] * x[j];
        opt = std::min(opt, obj);
      } else {
        const double v = max_lhs(p, x, negc, kInf);
        if (v == -kInf) return;
        opt = std::min(opt, -v);
      }
      pts.push_back(x);
    });
    if (pts.empty() || !is_finite(opt)) continue;
    const double cutoff = cutoff_on ? opt + 1e-7 : kInf;
    DomainPropagator prop(p);
    prop.continuous = nc > 0;
    if (cutoff_on) {
      prop.enable_objective();
      prop.set_objective_cutoff(cutoff);
    }
    BoundTrail trail;
    trail.init(n, 10000);
    ConflictPool pool(n, 200);
    auto valid = [&](const Conflict& c) {
      for (const auto& x : pts) {
        if (nc == 0) {
          double obj = 0.0;
          for (Index j = 0; j < n; ++j) obj += p.cost[j] * x[j];
          if (obj > cutoff) continue;
          if (c.violated_by(x, 1e-9)) return false;
        } else if (max_lhs(p, x, c, cutoff) > c.rhs + 1e-6 * (1.0 + std::abs(c.rhs))) {
          return false;
        }
      }
      return true;
    };
    for (int dive = 0; dive < 25; ++dive) {
      std::vector<double> lb = p.col_lb, ub = p.col_ub;
      trail.clear();
      prop.trail = &trail;
      bool ok = true;
      auto propagate = [&](std::vector<Index> seeds) {
        std::vector<Index> tight, ctight;
        for (int pass = 0; pass < 4 && ok && !seeds.empty(); ++pass) {
          tight.clear();
          ok = prop.propagate(lb, ub, seeds, tight);
          if (!ok) break;
          seeds.insert(seeds.end(), tight.begin(), tight.end());
          ctight.clear();
          ok = pool.propagate(p.is_int, p.col_lb, p.col_ub, lb, ub, seeds, ctight, dive, &trail);
          seeds = ctight;
        }
      };
      const bool bulk = dive % 2 == 1;
      std::vector<Index> decided;
      for (int level = 0; level < 8 && ok; ++level) {
        std::vector<Index> open;
        for (Index j = 0; j < n; ++j) {
          if (p.is_int[j] && lb[j] < ub[j]) open.push_back(j);
        }
        if (open.empty()) break;
        const Index j = open[rng() % open.size()];
        const double v = lb[j] + static_cast<int>(rng() % static_cast<unsigned>(ub[j] - lb[j] + 1));
        ++trail.level;
        if (rng() % 2) {
          if (v == ub[j]) continue;
          ub[j] = v;
          trail.push(j, true, v, BoundTrail::kDecision);
        } else {
          if (v == lb[j]) continue;
          lb[j] = v;
          trail.push(j, false, v, BoundTrail::kDecision);
        }
        if (bulk) decided.push_back(j);
        else propagate({j});
      }
      if (bulk && ok) propagate(decided);
      prop.trail = nullptr;
      if (ok) continue;
      ++st.calls;
      GenresOptions go;
      go.max_nnz = 50;
      GenresInfo info;
      Conflict g;
      const bool learned = analyze_genres(trail, prop.rows(), prop.row_lb(), prop.row_ub(), p.is_int, p.col_lb,
                                          p.col_ub, go, g, &info);
      st.steps += info.steps, st.reduced += info.reduced, st.mir += info.mir;
      std::vector<Conflict> out;
      analyze_propagation_conflict(trail, prop.rows(), prop.row_lb(), prop.row_ub(), p.is_int, p.col_lb, p.col_ub,
                                   50, 1, out);
      if (learned) {
        ++st.learned;
        st.nnz += static_cast<int>(g.idx.size());
        if (info.outcome == GenresInfo::kUip) ++st.uip;
        if (info.outcome == GenresInfo::kPartial) ++st.partial;
        if (!valid(g)) ++st.bad;
        // It must prove the node infeasible: minimal activity over the node's bounds above the rhs.
        double act = 0.0;
        for (size_t k = 0; k < g.idx.size(); ++k) act += g.val[k] * (g.val[k] > 0 ? lb[g.idx[k]] : ub[g.idx[k]]);
        if (!(act > g.rhs)) ++st.bad;
        pool.add(std::make_shared<Conflict>(std::move(g)), dive);
      } else {
        ++st.fails[std::min<int>(info.reason, 8)];
      }
      for (Conflict& c : out) pool.add(std::make_shared<Conflict>(std::move(c)), dive);
    }
  }
}

void print_genres(const char* what, const GenresStats& st) {
  std::printf("  %s: %d calls, %d learned (%d UIP, %d partial, %.1f nnz), %d steps, %d reduced, %d MIR, %d invalid;"
              " failures by reason",
              what, st.calls, st.learned, st.uip, st.partial,
              st.learned ? static_cast<double>(st.nnz) / st.learned : 0.0, st.steps, st.reduced, st.mir, st.bad);
  for (int f : st.fails) std::printf(" %d", f);
  std::printf("\n");
}

void set_env_value(const char* name, const char* value) {
#ifdef _WIN32
  _putenv_s(name, value ? value : "");
#else
  if (value) setenv(name, value, 1);
  else unsetenv(name);
#endif
}

}  // namespace

TEST(genres_conflicts_valid_pure_integer) {
  std::mt19937 rng(23);
  GenresStats st;
  genres_dives(rng, scaled(300), 0, false, st);
  print_genres("pure integer", st);
  CHECK(st.learned > 50);
  CHECK(st.mir > 0);
  CHECK_EQ(st.bad, 0);
  GenresStats sc;
  genres_dives(rng, scaled(200), 0, true, sc);
  print_genres("pure integer, cutoff at the optimum", sc);
  CHECK(sc.learned > 20);
  CHECK_EQ(sc.bad, 0);
}

TEST(genres_conflicts_valid_mixed) {
  std::mt19937 rng(29);
  GenresStats st;
  genres_dives(rng, scaled(120), 3, false, st);
  print_genres("mixed", st);
  CHECK(st.learned > 10);
  CHECK_EQ(st.bad, 0);
  GenresStats sc;
  genres_dives(rng, scaled(80), 3, true, sc);
  print_genres("mixed, cutoff at the optimum", sc);
  CHECK_EQ(sc.bad, 0);
}

TEST(genres_hand_example) {
  // x0 + x1 + x2 >= 2 (row 0), x0 + x1 + x3 <= 1 (row 1), binaries. The decision x2 <= 0 implies
  // x0 >= 1 and x1 >= 1 by row 0, and row 1 fails. Resolution eliminates x1 and x0; the learned
  // constraint must hold for every feasible point and fail under x2 <= 0.
  ModelBuilder b;
  for (int j = 0; j < 4; ++j) b.add_col("x" + std::to_string(j), 0.0, 0.0, 1.0, VarType::Binary);
  const Index r0 = b.add_row("r0", 2.0, kInf);
  for (int j : {0, 1, 2}) b.add_coef(r0, j, 1.0);
  const Index r1 = b.add_row("r1", -kInf, 1.0);
  for (int j : {0, 1, 3}) b.add_coef(r1, j, 1.0);
  MipProblem p;
  REQUIRE(MipProblem::from_model(b.build(), &p));
  DomainPropagator prop(p);
  BoundTrail trail;
  trail.init(4, 100);
  prop.trail = &trail;
  std::vector<double> lb = p.col_lb, ub = p.col_ub;
  ++trail.level;
  ub[2] = 0.0;
  trail.push(2, true, 0.0, BoundTrail::kDecision);
  std::vector<Index> t;
  CHECK(!prop.propagate(lb, ub, {2}, t));
  GenresInfo info;
  Conflict g;
  REQUIRE(analyze_genres(trail, prop.rows(), prop.row_lb(), prop.row_ub(), p.is_int, p.col_lb, p.col_ub,
                         GenresOptions(), g, &info));
  CHECK(info.steps >= 1);
  int feasible = 0;
  for (int mask = 0; mask < 16; ++mask) {
    std::vector<double> x(4);
    for (int j = 0; j < 4; ++j) x[j] = (mask >> j) & 1;
    if (x[0] + x[1] + x[2] < 2 || x[0] + x[1] + x[3] > 1) continue;
    ++feasible;
    CHECK(!g.violated_by(x, 1e-9));
  }
  CHECK(feasible > 0);
  std::vector<double> l2 = p.col_lb, u2 = p.col_ub;
  u2[2] = 0.0;
  CHECK(!propagate_conflict(g, p.is_int, l2, u2, nullptr, nullptr));
}

TEST(genres_end_to_end_same_optimum) {
  // Cut-based conflicts (with the clauses, and instead of them) must not change the optimum.
  std::mt19937 rng(78);
  SolverParams params;
  params.verbosity = 0;
  params.mip_rel_gap = 1e-9;
  params.mip_abs_gap = 1e-9;
  params.node_limit = 20000;
  params.presolve = false;
  int compared = 0;
  for (int t = 0; t < scaled(30); ++t) {
    Model m = random_mip(rng, 25 + rng() % 15, rng() % 6, 10 + rng() % 8);
    params.threads = 1 + t % 2;
    params.deterministic = t % 4 == 1;
    set_env("MILP_NO_CONFLICT", true);
    SolveResult r0 = solve(m, params);
    set_env("MILP_NO_CONFLICT", false);
    set_env_value("MILP_GENRES", "1");
    SolveResult r1 = solve(m, params);
    set_env_value("MILP_GENRES", "2");
    SolveResult r2 = solve(m, params);
    set_env_value("MILP_GENRES", nullptr);
    if (r0.status == Status::Infeasible) {
      CHECK(r1.status == Status::Infeasible);
      CHECK(r2.status == Status::Infeasible);
      continue;
    }
    if (r0.status != Status::Optimal || r1.status != Status::Optimal || r2.status != Status::Optimal) continue;
    CHECK_NEAR(r1.objective, r0.objective, 1e-6 * std::max(1.0, std::abs(r0.objective)));
    CHECK_NEAR(r2.objective, r0.objective, 1e-6 * std::max(1.0, std::abs(r0.objective)));
    ++compared;
  }
  std::printf("  %d instances compared\n", compared);
  CHECK(compared > 10);
}
