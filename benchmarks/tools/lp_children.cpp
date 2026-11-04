// Child-LP experiment (v0.6 node-LP design, PLAN.md 2026-11-03): root LP of an MPS model (integrality
// relaxed), then warm-started re-solves after one bound change (both branches of the k most fractional
// columns, each from the saved root state) and a dive (always the most fractional column, warm from the
// previous LP). The candidates file lets lp_children_highs.py replay the same bound changes in HiGHS.
// Not part of the build:
//   c++ -std=c++20 -O2 -Isrc benchmarks/tools/lp_children.cpp build/libmilp_lib.a -o /tmp/lp_children
//   MILP_PRESOLVE_ONLY=1 MILP_PRESOLVE_DUMP=/tmp/m.mps build/popso <model.mps>     # our presolved model
//   /tmp/lp_children /tmp/m.mps 10 10 /tmp/m.cands      # VERB=3: per-100-iteration shifts / zero steps / flips
//   .venv/bin/python benchmarks/tools/lp_children_highs.py /tmp/m.mps /tmp/m.cands
// Usage: lp_children <model.mps> <k candidates> <dive depth> [candidates file out]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "milp/io/mps_reader.hpp"
#include "milp/lp/simplex.hpp"
#include "milp/search/mip_problem.hpp"

using namespace milp;

static double now() {
  using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();
}

int main(int argc, char** argv) {
  Model model = read_mps(argv[1]);
  MipProblem p;
  if (!MipProblem::from_model(model, &p)) return 1;
  const int K = argc > 2 ? std::atoi(argv[2]) : 20;
  const int depth = argc > 3 ? std::atoi(argv[3]) : 20;
  FILE* cf = argc > 4 ? std::fopen(argv[4], "w") : nullptr;
  DualSimplex lp;
  lp.load(p.A, p.cost, p.col_lb, p.col_ub, p.row_lb, p.row_ub);
  lp.set_seed(0);
  LpOptions opt;
  opt.verbosity = std::getenv("VERB") ? std::atoi(std::getenv("VERB")) : 0;
  double t = now();
  LpStatus st = lp.solve(opt);
  std::printf("root: %s obj %.10g its %lld time %.2f (m %d n %d)\n", to_string(st), lp.objective() + p.offset,
              lp.iterations(), now() - t, p.num_rows(), p.num_cols());
  DualSimplex::State root;
  lp.save_state(root);
  std::vector<double> x = lp.col_values();
  std::vector<std::pair<double, Index>> frac;
  for (Index j = 0; j < p.num_cols(); ++j) {
    if (!p.is_int[j]) continue;
    const double f = x[j] - std::floor(x[j]);
    if (f > 1e-6 && f < 1 - 1e-6) frac.push_back({-std::min(f, 1 - f), j});
  }
  std::sort(frac.begin(), frac.end());
  std::printf("fractional: %zu\n", frac.size());
  long long tot = 0;
  double ttot = 0;
  int cnt = 0;
  for (int k = 0; k < K && k < (int)frac.size(); ++k) {
    const Index j = frac[k].second;
    for (int dir = 0; dir < 2; ++dir) {
      lp.restore_state(root);
      const double lb = p.col_lb[j], ub = p.col_ub[j];
      if (dir == 0) lp.set_col_bounds(j, lb, std::floor(x[j]));
      else lp.set_col_bounds(j, std::ceil(x[j]), ub);
      const long long it0 = lp.iterations();
      t = now();
      st = lp.solve(opt);
      const double dt = now() - t;
      std::printf("child col %d %s: %s obj %.10g its %lld time %.3f\n", j, dir ? "up" : "down", to_string(st),
                  lp.objective() + p.offset, lp.iterations() - it0, dt);
      if (cf) std::fprintf(cf, "%d %d %.17g\n", j, dir, dir ? std::ceil(x[j]) : std::floor(x[j]));
      tot += lp.iterations() - it0;
      ttot += dt;
      ++cnt;
      lp.set_col_bounds(j, lb, ub);
    }
  }
  if (cnt) std::printf("children: %d, mean its %.0f, mean time %.3f\n", cnt, double(tot) / cnt, ttot / cnt);
  // Dive: always branch down/up (alternating by fraction) on the most fractional column, warm from the previous LP.
  lp.restore_state(root);
  std::vector<double> lbv = p.col_lb, ubv = p.col_ub;
  for (int d = 0; d < depth; ++d) {
    std::vector<double> xc = lp.col_values();
    Index best = -1;
    double bf = 0;
    for (Index j = 0; j < p.num_cols(); ++j) {
      if (!p.is_int[j]) continue;
      const double f = xc[j] - std::floor(xc[j]);
      const double s = std::min(f, 1 - f);
      if (s > 1e-6 && s > bf) bf = s, best = j;
    }
    if (best < 0) break;
    const double f = xc[best] - std::floor(xc[best]);
    if (f < 0.5) ubv[best] = std::floor(xc[best]);
    else lbv[best] = std::ceil(xc[best]);
    lp.set_col_bounds(best, lbv[best], ubv[best]);
    const long long it0 = lp.iterations();
    t = now();
    st = lp.solve(opt);
    std::printf("dive %d col %d %s: %s obj %.10g its %lld time %.3f\n", d, best, f < 0.5 ? "down" : "up",
                to_string(st), lp.objective() + p.offset, lp.iterations() - it0, now() - t);
    if (cf) std::fprintf(cf, "D %d %.17g %.17g\n", best, lbv[best], ubv[best]);
    if (st != LpStatus::Optimal) break;
  }
  if (cf) std::fclose(cf);
  return 0;
}
