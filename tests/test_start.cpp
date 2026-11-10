// Tests for user starting solutions (solver/start.hpp, solve(model, params, start)).

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "milp/core/log.hpp"
#include "milp/core/timer.hpp"
#include "milp/io/mps_reader.hpp"
#include "milp/search/branch_and_bound.hpp"
#include "milp/search/mip_problem.hpp"
#include "milp/solver/solver.hpp"
#include "milp/solver/start.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// max 5x + 4y + 3z  s.t.  2x + 3y + z <= 5,  4x + y + 2z <= 11,  3x + 4y + 2z <= 8,  x, y, z >= 0 integer,
// plus a continuous slack-like column w <= x + y with objective 0.5 (w in [0, 10]). Optimum 14 at (2, 0, 1, 2).
Model small_mip() {
  ModelBuilder b;
  b.model().sense = ObjSense::Maximize;
  const Index x = b.add_col("x", 5, 0, kInf, VarType::Integer);
  const Index y = b.add_col("y", 4, 0, kInf, VarType::Integer);
  const Index z = b.add_col("z", 3, 0, kInf, VarType::Integer);
  const Index w = b.add_col("w", 0.5, 0, 10);
  const double a[3][3] = {{2, 3, 1}, {4, 1, 2}, {3, 4, 2}};
  const double rhs[3] = {5, 11, 8};
  for (int i = 0; i < 3; ++i) {
    const Index r = b.add_row("r" + std::to_string(i), -kInf, rhs[i]);
    b.add_coef(r, x, a[i][0]);
    b.add_coef(r, y, a[i][1]);
    b.add_coef(r, z, a[i][2]);
  }
  const Index r = b.add_row("link", -kInf, 0);
  b.add_coef(r, w, 1);
  b.add_coef(r, x, -1);
  b.add_coef(r, y, -1);
  return b.build();
}

SolverParams quiet_params() {
  SolverParams p;
  p.verbosity = 0;
  p.threads = 1;
  return p;
}

StartOutcome prepare(const Model& m, const std::vector<double>& start) {
  MipProblem p;
  REQUIRE(MipProblem::from_model(m, &p));
  Timer t;
  return prepare_start(p, start, quiet_params(), t);
}

}  // namespace

TEST(start_complete_feasible_is_kept) {
  const Model m = small_mip();
  const std::vector<double> s = {1, 1, 0, 2};  // objective 5 + 4 + 1 = 10
  const StartOutcome o = prepare(m, s);
  REQUIRE(o.x.size() == 4);
  CHECK(o.x == s);
  CHECK(o.message.find("feasible") != std::string::npos);
}

TEST(start_complete_infeasible_is_repaired) {
  const Model m = small_mip();
  // Integers fine, continuous w violates the link row: the LP over w repairs it.
  const StartOutcome o = prepare(m, {1, 1, 0, 7});
  REQUIRE(o.x.size() == 4);
  CHECK_NEAR(o.x[0], 1, 1e-9);
  CHECK_NEAR(o.x[1], 1, 1e-9);
  CHECK(m.violation(o.x).max() <= 1e-6);
}

TEST(start_partial_is_completed) {
  const Model m = small_mip();
  const StartOutcome o = prepare(m, {2, kNaN, kNaN, kNaN});
  REQUIRE(o.x.size() == 4);
  CHECK_NEAR(o.x[0], 2, 1e-9);  // the given integer stays fixed
  CHECK(m.violation(o.x).max() <= 1e-6);
  CHECK(o.message.find("completed") != std::string::npos);
}

TEST(start_rejected_cases) {
  const Model m = small_mip();
  // x = 2, y = 2 violates 2x + 3y + z <= 5 whatever z, w are.
  StartOutcome o = prepare(m, {2, 2, kNaN, kNaN});
  CHECK(o.x.empty());
  CHECK(o.message.find("ignored") != std::string::npos);
  o = prepare(m, {kNaN, kNaN, kNaN, 1.0});  // no integer values
  CHECK(o.x.empty());
  o = prepare(m, {kNaN, kNaN, kNaN, kNaN});
  CHECK(o.x.empty());
  o = prepare(m, {1, 1});  // wrong size
  CHECK(o.x.empty());
}

TEST(start_solve_end_to_end) {
  const Model m = small_mip();
  for (int variant = 0; variant < 4; ++variant) {
    SolverParams p = quiet_params();
    p.presolve = variant % 2 == 0;
    p.deterministic = variant >= 2;
    p.threads = variant >= 2 ? 2 : 1;
    for (const std::vector<double>& s : {std::vector<double>{1, 1, 0, 2}, std::vector<double>{2, kNaN, kNaN, kNaN},
                                         std::vector<double>{2, 2, kNaN, kNaN}, std::vector<double>{2, 0, 1, 2}}) {
      const SolveResult r = solve(m, p, s);
      CHECK(r.status == Status::Optimal);
      REQUIRE(r.has_solution);
      CHECK_NEAR(r.objective, 14.0, 1e-6);
      CHECK_NEAR(r.bound, 14.0, 1e-6);
    }
  }
}

TEST(start_is_incumbent_under_node_limit) {
  // The known optimum of a benchmark instance as the start, with a search that stops after
  // one node: the result keeps (at least) the start's objective.
  const Model m = read_mps(std::string(MILP_SOURCE_DIR) + "/benchmarks/instances/knapsack_small_1.mps");
  SolverParams p = quiet_params();
  const SolveResult ref = solve(m, p);
  REQUIRE(ref.status == Status::Optimal);
  p.node_limit = 1;
  const SolveResult r = solve(m, p, ref.x);
  REQUIRE(r.has_solution);
  const double sense = m.sense == ObjSense::Maximize ? -1.0 : 1.0;
  CHECK(sense * r.objective <= sense * ref.objective + 1e-6 * (1.0 + std::abs(ref.objective)));
  CHECK(m.violation(r.x).max() <= 1e-6);
  // Half of the integer values given: completed (or ignored), never an invalid result.
  std::vector<double> half(m.num_cols(), kNaN);
  for (Index j = 0; j < m.num_cols(); j += 2) half[j] = ref.x[j];
  const SolveResult h = solve(m, p, half);
  if (h.has_solution) CHECK(m.violation(h.x).max() <= 1e-6);
}

TEST(start_rejected_cutoff_from_its_own_objective) {
  // A restart hands the search the incumbent before it, crushed through the restart's presolve, and its
  // objective (SearchOptions::start_objective). A dual reduction may remove that point: the crushed point
  // is then rejected, and the cutoff must come from the start's objective, not from the crushed point's
  // (here 1, below the optimum 3: nothing would be found, and the caller would return its old incumbent
  // as optimal; neos-3083819-nubu on x86).
  // min x + y  s.t.  x + y >= 3,  x, y integer in [0, 10]: optimum 3.
  MipProblem p;
  {
    ModelBuilder b;
    const Index x = b.add_col("x", 1, 0, 10, VarType::Integer);
    const Index y = b.add_col("y", 1, 0, 10, VarType::Integer);
    const Index r = b.add_row("r", 3, kInf);
    b.add_coef(r, x, 1);
    b.add_coef(r, y, 1);
    REQUIRE(MipProblem::from_model(b.build(), &p));
  }
  SolverParams sp = quiet_params();
  Logger log(0);
  for (const bool with_objective : {true, false}) {
    SearchOptions opts;
    opts.check_resident = false;
    opts.feasibility_jump = false;
    opts.start = {1.0, 0.0};  // violates the row; its objective here is 1
    if (with_objective) opts.start_objective = 5.0;  // the incumbent's own objective
    Timer timer;
    const MipResult r = branch_and_bound(p, sp, log, timer, opts);
    if (with_objective) {
      CHECK(r.status == Status::Optimal);
      REQUIRE(r.has_solution);
      CHECK_NEAR(r.objective, 3.0, 1e-9);
    } else {
      CHECK(!r.has_solution);  // the point's own objective as the cutoff: below the optimum
    }
  }
}

TEST(start_repaired_after_presolve) {
  // A start that violates the (presolved) problem in one row: the repair keeps its integers outside that row (z, w)
  // and completes the rest by a sub-MIP, which becomes the first incumbent ("start repair"), in both modes.
  // min x + y + z + 2 w  s.t.  x + y >= 3,  z + w >= 1,  integers in [0, 10]: optimum 4 (z = 1, w = 0).
  MipProblem p;
  {
    ModelBuilder b;
    const Index x = b.add_col("x", 1, 0, 10, VarType::Integer);
    const Index y = b.add_col("y", 1, 0, 10, VarType::Integer);
    const Index z = b.add_col("z", 1, 0, 10, VarType::Integer);
    const Index w = b.add_col("w", 2, 0, 10, VarType::Integer);
    const Index r = b.add_row("r", 3, kInf);
    b.add_coef(r, x, 1);
    b.add_coef(r, y, 1);
    const Index s = b.add_row("s", 1, kInf);
    b.add_coef(s, z, 1);
    b.add_coef(s, w, 1);
    REQUIRE(MipProblem::from_model(b.build(), &p));
  }
  for (const bool det : {false, true}) {
    SolverParams sp = quiet_params();
    sp.deterministic = det;
    Logger log(0);
    SearchOptions opts;
    opts.check_resident = false;
    opts.feasibility_jump = false;
    opts.start = {1.0, 0.0, 1.0, 0.0};  // violates r
    opts.start_objective = 9.0;         // the user's start was feasible before presolve, at 9
    std::vector<std::string> sources;
    opts.on_trace = [&](double, const char* source, double) { sources.emplace_back(source); };
    Timer timer;
    const MipResult r = branch_and_bound(p, sp, log, timer, opts);
    CHECK(r.status == Status::Optimal);
    REQUIRE(r.has_solution);
    CHECK_NEAR(r.objective, 4.0, 1e-9);
    REQUIRE(!sources.empty());
    CHECK(sources.front() == "start repair");
  }
}
