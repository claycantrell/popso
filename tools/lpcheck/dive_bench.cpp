// LP engine benchmark / warm-start checker on a MIP: simulates the node LP
// sequence of a depth-first search (bound changes + warm solves, backtracking
// via restore_state or set_basis, occasional strong-branching style probes).
//   ./dive_bench <model.mps> [solves] [check_every]
// check_every > 0 re-solves every k-th LP from scratch and compares status and
// objective (a mismatch is printed and counted).
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>

#include "milp/io/mps_reader.hpp"
#include "milp/lp/simplex.hpp"

using namespace milp;

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: dive_bench <model.mps> [solves] [check_every]\n");
    return 1;
  }
  Model mdl = read_mps(argv[1]);
  const long long solves = argc > 2 ? std::atoll(argv[2]) : 2000;
  const long long check_every = argc > 3 ? std::atoll(argv[3]) : 0;
  const Index n = mdl.num_cols();
  std::vector<double> cost = mdl.obj;
  if (mdl.sense == ObjSense::Maximize)
    for (double& c : cost) c = -c;
  DualSimplex lp;
  lp.load(mdl.A, cost, mdl.col_lb, mdl.col_ub, mdl.row_lb, mdl.row_ub);
  LpOptions opt;
  auto t0 = std::chrono::steady_clock::now();
  LpStatus st = lp.solve(opt);
  const double root_obj = lp.objective();
  const double t_root = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  std::printf("root: %s obj %.10g, %lld iterations, %.3fs\n", to_string(st), root_obj, lp.iterations(), t_root);
  if (st != LpStatus::Optimal) return 0;

  struct Level {
    Index j;
    double lb, ub;  // bounds before the branch
    DualSimplex::State state;
  };
  std::vector<Level> stack;
  std::mt19937 rng(7);
  long long done = 0, mismatches = 0, probes = 0;
  const long long it0 = lp.iterations();
  t0 = std::chrono::steady_clock::now();
  DualSimplex::State probe_state;
  while (done < solves) {
    // Pick a random fractional integer column.
    Index best = -1;
    if (st == LpStatus::Optimal) {
      const double* x = lp.x();
      long long count = 0;
      for (Index j = 0; j < n; ++j) {
        if (mdl.col_type[j] == VarType::Continuous || std::abs(x[j] - std::round(x[j])) <= 1e-6) continue;
        if (rng() % ++count == 0) best = j;
      }
    }
    if (best < 0) {
      // Backtrack a random number of levels (all of them half of the time, so
      // the dives restart from the root).
      if (stack.empty()) break;
      const size_t k = rng() % 2 ? stack.size() : 1 + rng() % std::min<size_t>(stack.size(), 8);
      for (size_t t = 0; t < k; ++t) {
        Level& L = stack.back();
        lp.set_col_bounds(L.j, L.lb, L.ub);
        if (t + 1 == k) {
          if (rng() % 3 == 0) lp.set_basis(Basis{L.state.status});
          else lp.restore_state(L.state);
        }
        stack.pop_back();
      }
    } else {
      const double v = lp.x()[best];
      if (rng() % 8 == 0) {
        // Strong-branching style probe with an iteration limit.
        lp.save_state(probe_state);
        const double ol = lp.col_lower(best), ou = lp.col_upper(best);
        lp.set_col_bounds(best, std::ceil(v), ou);
        LpOptions po = opt;
        po.iteration_limit = 30;
        lp.solve(po);
        lp.set_col_bounds(best, ol, ou);
        lp.restore_state(probe_state);
        ++probes;
      }
      Level L{best, lp.col_lower(best), lp.col_upper(best), {}};
      lp.save_state(L.state);
      stack.push_back(std::move(L));
      if (rng() % 2) lp.set_col_bounds(best, std::ceil(v), stack.back().ub);
      else lp.set_col_bounds(best, stack.back().lb, std::floor(v));
    }
    st = lp.solve(opt);
    ++done;
    if (check_every > 0 && done % check_every == 0) {
      std::vector<double> lb(n), ub(n);
      for (Index j = 0; j < n; ++j) lb[j] = lp.col_lower(j), ub[j] = lp.col_upper(j);
      DualSimplex fresh;
      fresh.load(mdl.A, cost, lb, ub, mdl.row_lb, mdl.row_ub);
      const LpStatus st2 = fresh.solve(opt);
      bool bad = st != st2;
      if (!bad && st == LpStatus::Optimal)
        bad = std::abs(fresh.objective() - lp.objective()) > 1e-6 * (1 + std::abs(lp.objective()));
      if (bad) {
        ++mismatches;
        std::printf("MISMATCH solve %lld: warm %s %.10g fresh %s %.10g\n", done, to_string(st), lp.objective(),
                    to_string(st2), fresh.objective());
      }
    }
  }
  const double t = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  const long long its = lp.iterations() - it0;
  std::printf("solves %lld (+%lld probes), iterations %lld, time %.3fs, %.1f us/solve, %.2f us/iter, mismatches %lld\n",
              done, probes, its, t, 1e6 * t / done, its ? 1e6 * t / its : 0.0, mismatches);
}
