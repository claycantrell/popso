// Random LP KKT checker for DualSimplex.
#include "milp/lp/simplex.hpp"
#include <cstdio>
#include <random>
using namespace milp;
int main(int argc, char** argv) {
  int trials = argc > 1 ? atoi(argv[1]) : 500;
  int big = argc > 2 ? atoi(argv[2]) : 1;  // size multiplier: larger LPs exercise the LU kernel / dense block
  int dens = big == 1 ? 25 : 25 / big + 2;
  std::mt19937 rng(7);
  int counts[8] = {0}; int bad = 0;
  for (int t = 0; t < trials; ++t) {
    int m = 1 + rng() % (40 * big), n = 1 + rng() % (60 * big);
    std::vector<Index> mj, mi; std::vector<double> mv;
    for (int j = 0; j < n; ++j) for (int i = 0; i < m; ++i) if ((int)(rng() % 100) < dens) { mj.push_back(j); mi.push_back(i); mv.push_back((int)(rng() % 21) - 10); }
    auto A = SparseMatrix::from_triplets(n, m, mj, mi, mv);
    std::vector<double> c(n), cl(n), cu(n), rl(m), ru(m);
    auto pick = [&](double& l, double& u, int scale) {
      int k = rng() % 6; double a = (int)(rng() % (2*scale+1)) - scale, b = a + rng() % (scale+1);
      if (k == 0) { l = a; u = b; } else if (k == 1) { l = a; u = kInf; } else if (k == 2) { l = -kInf; u = b; }
      else if (k == 3) { l = -kInf; u = kInf; } else if (k == 4) { l = u = a; } else { l = 0; u = kInf; }
    };
    for (int j = 0; j < n; ++j) { c[j] = (int)(rng() % 21) - 10; pick(cl[j], cu[j], 5); if (t % 3 == 0) { if (!std::isfinite(cl[j])) cl[j] = -10; if (!std::isfinite(cu[j])) cu[j] = 10; } }
    for (int i = 0; i < m; ++i) pick(rl[i], ru[i], 20);
    DualSimplex lp; lp.load(A, c, cl, cu, rl, ru);
    LpOptions opt; opt.perturb = (t % 2 == 0);
    LpStatus st = lp.solve(opt);
    counts[(int)st]++;
    { FILE* f = fopen("lps.txt", t == 0 ? "w" : "a"); fprintf(f, "LP %d %d %d %.17g\n", n, m, (int)st, st == LpStatus::Optimal ? lp.objective() : 0.0);
      for (int j = 0; j < n; ++j) fprintf(f, "%g %g %g\n", c[j], cl[j], cu[j]);
      for (int i = 0; i < m; ++i) fprintf(f, "%g %g\n", rl[i], ru[i]);
      fprintf(f, "%d\n", A.nnz()); for (int j = 0; j < n; ++j) for (Index p = A.begin(j); p < A.end(j); ++p) fprintf(f, "%d %d %g\n", j, A.index[p], A.value[p]);
      fclose(f); }
    if (st == LpStatus::Optimal) {
      auto x = lp.col_values(); auto act = lp.row_activities(); auto y = lp.row_duals(); auto d = lp.reduced_costs();
      double pinf = 0, dinf = 0, err = 0;
      std::vector<double> ax(m, 0);
      for (int j = 0; j < n; ++j) for (Index p = A.begin(j); p < A.end(j); ++p) ax[A.index[p]] += A.value[p] * x[j];
      for (int i = 0; i < m; ++i) { err = std::max(err, std::abs(ax[i] - act[i])); pinf = std::max({pinf, rl[i] - ax[i], ax[i] - ru[i]}); }
      for (int j = 0; j < n; ++j) pinf = std::max({pinf, cl[j] - x[j], x[j] - cu[j]});
      // dual feasibility / complementary slackness for columns
      for (int j = 0; j < n; ++j) {
        bool atl = std::abs(x[j] - cl[j]) < 1e-7, atu = std::abs(x[j] - cu[j]) < 1e-7;
        if (d[j] > 1e-6 && !atl) dinf = std::max(dinf, d[j]);
        if (d[j] < -1e-6 && !atu) dinf = std::max(dinf, -d[j]);
      }
      // rows: logical s=-ax with d_s = -y
      for (int i = 0; i < m; ++i) {
        double ds = -y[i]; double s = -ax[i];
        bool atl = std::abs(s - (-ru[i])) < 1e-7, atu = std::abs(s - (-rl[i])) < 1e-7;
        if (ds > 1e-6 && !atl) dinf = std::max(dinf, ds);
        if (ds < -1e-6 && !atu) dinf = std::max(dinf, -ds);
      }
      if (pinf > 1e-6 || dinf > 1e-6 || err > 1e-6) { ++bad; printf("trial %d m=%d n=%d pinf %g dinf %g err %g\n", t, m, n, pinf, dinf, err); }
    }
  }
  printf("optimal %d infeasible %d unbounded %d itlim %d numerr %d | KKT failures %d\n", counts[0], counts[1], counts[2], counts[4], counts[6], bad);
}
