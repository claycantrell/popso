#include "milp/lp/lu.hpp"
#include <cstdio>
#include <cstdlib>
#include <random>
#include <cmath>
using namespace milp;
static int g_updates = 30;

// Large sparse bases (LP-like: many unit columns, short structural columns, optional dependent
// columns), with eta updates; checks dense and HVector solves by sparse residuals.
static double large_trials(std::mt19937& rng, int trials) {
  double worst = 0;
  for (int trial = 0; trial < trials; ++trial) {
    int m = 50 + rng() % 1500;
    int slack_pct = rng() % 90, len = 1 + rng() % 6;
    std::vector<std::vector<std::pair<Index, double>>> cols(m);
    for (int c = 0; c < m; ++c) {
      if ((int)(rng() % 100) < slack_pct) { cols[c].push_back({(Index)(rng() % m), 1.0}); continue; }
      std::vector<char> used(m, 0);
      for (int k = 0; k < len; ++k) {
        Index r = rng() % m;
        if (used[r]) continue;
        used[r] = 1;
        cols[c].push_back({r, (double)((int)(rng() % 19) - 9) + ((rng() % 4) ? 0.0 : 0.37)});
      }
    }
    if (trial % 5 == 0 && m > 3) cols[2] = cols[1];  // dependent column
    if (trial % 7 == 0 && m > 3) cols[3].clear();    // empty column
    auto build = [&](std::vector<Index>& st, std::vector<Index>& ix, std::vector<double>& vv) {
      st.assign(1, 0); ix.clear(); vv.clear();
      for (int c = 0; c < m; ++c) { for (auto [r, v] : cols[c]) { ix.push_back(r); vv.push_back(v); } st.push_back(ix.size()); }
    };
    std::vector<Index> st, ix; std::vector<double> vv;
    build(st, ix, vv);
    LuFactor lu; std::vector<std::pair<Index, Index>> rep;
    if (getenv("THR")) lu.pivot_threshold = atof(getenv("THR"));
    lu.factorize(m, st, ix, vv, &rep);
    for (auto [p, r] : rep) cols[p] = {{r, 1.0}};
    for (int u = 0; u < g_updates; ++u) {
      std::vector<std::pair<Index, double>> a;
      for (int k = 0; k < 1 + (int)(rng() % 5); ++k) a.push_back({(Index)(rng() % m), (double)((int)(rng() % 7) - 3)});
      std::vector<double> alpha(m, 0); for (auto [r, v] : a) alpha[r] += v;
      std::vector<std::pair<Index, double>> ac; for (int r = 0; r < m; ++r) if (alpha[r] != 0) ac.push_back({r, alpha[r]});
      HVector h; h.resize(m); for (auto [r, v] : ac) { h.dense[r] = v; h.nz[h.count++] = r; }
      lu.ftran(h);
      int p = rng() % m;
      double amax = 0; for (int i = 0; i < m; ++i) amax = std::max(amax, std::abs(h.dense[i]));
      if (std::abs(h.dense[p]) < std::max(0.1 * amax, 1e-3)) continue;  // simplex-like (well conditioned) pivots only
      if (u % 2) lu.update(p, h); else lu.update(p, h.dense);
      cols[p] = ac;
    }
    if (getenv("FRESH")) {  // refactor the updated basis (separates LU from eta-update error)
      build(st, ix, vv);
      lu.factorize(m, st, ix, vv, &rep);
      for (auto [p, r] : rep) cols[p] = {{r, 1.0}};
    }
    auto mulB = [&](const std::vector<double>& x) { std::vector<double> b(m, 0); for (int c = 0; c < m; ++c) for (auto [r, v] : cols[c]) b[r] += v * x[c]; return b; };
    auto mulBt = [&](const std::vector<double>& x) { std::vector<double> b(m, 0); for (int c = 0; c < m; ++c) for (auto [r, v] : cols[c]) b[c] += v * x[r]; return b; };
    for (int s = 0; s < 4; ++s) {
      // Sparse or dense right-hand sides.
      std::vector<double> rhs(m, 0);
      if (s < 2) { for (int k = 0; k < 3; ++k) rhs[rng() % m] = (int)(rng() % 9) - 4; }
      else for (int i = 0; i < m; ++i) rhs[i] = (int)(rng() % 11) - 5;
      std::vector<double> y = rhs; HVector h; h.resize(m); h.dense = rhs; h.rebuild_index();
      if (s % 2 == 0) { lu.ftran(y); lu.ftran(h); } else { lu.btran(y); lu.btran(h); }
      auto back = s % 2 == 0 ? mulB(y) : mulBt(y);
      double scale = 1, e = 0, eh = 0;
      for (int i = 0; i < m; ++i) scale = std::max(scale, std::abs(y[i]));
      for (int i = 0; i < m; ++i) e = std::max(e, std::abs(back[i] - rhs[i]));
      std::vector<char> in(m, 0);
      for (Index k = 0; k < h.count; ++k) in[h.nz[k]] = 1;
      for (int i = 0; i < m; ++i) { eh = std::max(eh, std::abs(h.dense[i] - y[i])); if (!in[i] && h.dense[i] != 0) eh = 1e300; }
      e /= scale;
      worst = std::max({worst, e, eh});
      if (e > 1e-9 || eh > 1e-9) printf("large trial %d m=%d rep=%zu %s rel err %g hvec err %g\n", trial, m, rep.size(), s % 2 ? "btran" : "ftran", e, eh);
    }
    // Multi-rhs ftran must agree with single ftrans.
    HVector h[3], g[3];
    for (int k = 0; k < 3; ++k) {
      h[k].resize(m);
      for (int t = 0; t < (k == 2 ? m : 3); ++t) h[k].dense[rng() % m] = (int)(rng() % 9) - 4;
      h[k].rebuild_index();
      g[k] = h[k];
      lu.ftran(g[k]);
    }
    HVector h2[2] = {h[0], h[1]};
    lu.ftran(h[0], h[1], h[2]);
    lu.ftran(h2[0], h2[1]);
    double em = 0;
    for (int k = 0; k < 3; ++k) for (int i = 0; i < m; ++i) em = std::max(em, std::abs(h[k].dense[i] - g[k].dense[i]));
    for (int k = 0; k < 2; ++k) for (int i = 0; i < m; ++i) em = std::max(em, std::abs(h2[k].dense[i] - g[k].dense[i]));
    if (em > 1e-12 * std::max(1.0, em)) printf("large trial %d multi-rhs ftran mismatch %g\n", trial, em);
    worst = std::max(worst, em > 1e-9 ? em : 0.0);
  }
  return worst;
}

int main(int argc, char** argv) {
  if (argc > 1) g_updates = atoi(argv[1]);
  std::mt19937 rng(1);
  double worst = 0;
  for (int trial = 0; trial < 300; ++trial) {
    int m = 1 + rng() % 60;
    double dens = 0.02 + (rng() % 100) / 300.0;
    std::vector<std::vector<double>> B(m, std::vector<double>(m, 0));
    for (int c = 0; c < m; ++c) {
      if (rng() % 3 == 0) { B[c][c] = 1; continue; } // slack-like column
      for (int r = 0; r < m; ++r) if ((rng() % 1000) / 1000.0 < dens) B[r][c] = (int)(rng() % 19) - 9;
      B[rng() % m][c] += 1 + rng() % 5;
    }
    if (trial % 7 == 0 && m > 2) for (int r = 0; r < m; ++r) B[r][1] = 2 * B[r][0]; // force singular
    std::vector<Index> st{0}, ix; std::vector<double> vv;
    for (int c = 0; c < m; ++c) { for (int r = 0; r < m; ++r) if (B[r][c] != 0) { ix.push_back(r); vv.push_back(B[r][c]); } st.push_back(ix.size()); }
    LuFactor lu; std::vector<std::pair<Index,Index>> rep;
    lu.factorize(m, st, ix, vv, &rep);
    for (auto [p, r] : rep) { for (int i = 0; i < m; ++i) B[i][p] = 0; B[r][p] = 1; }
    // updates: replace a few columns with random ones
    for (int u = 0; u < 5; ++u) {
      std::vector<double> a(m, 0); for (int r = 0; r < m; ++r) if (rng() % 3 == 0) a[r] = (int)(rng() % 7) - 3;
      std::vector<double> alpha = a; lu.ftran(alpha);
      int p = rng() % m; if (std::abs(alpha[p]) < 1e-3) continue;
      lu.update(p, alpha); for (int r = 0; r < m; ++r) B[r][p] = a[r];
    }
    std::vector<double> x(m), b(m, 0);
    for (int i = 0; i < m; ++i) x[i] = (int)(rng() % 11) - 5;
    for (int r = 0; r < m; ++r) for (int c = 0; c < m; ++c) b[r] += B[r][c] * x[c];
    lu.ftran(b); double e1 = 0; for (int i = 0; i < m; ++i) e1 = std::max(e1, std::abs(b[i] - x[i]));
    std::vector<double> e(m, 0); for (int c = 0; c < m; ++c) for (int r = 0; r < m; ++r) e[c] += B[r][c] * x[r];
    lu.btran(e); double e2 = 0; for (int i = 0; i < m; ++i) e2 = std::max(e2, std::abs(e[i] - x[i]));
    worst = std::max({worst, e1, e2});
    if (e1 > 1e-6 || e2 > 1e-6) printf("trial %d m=%d rep=%zu ftran err %g btran err %g\n", trial, m, rep.size(), e1, e2);
  }
  printf("worst error %g\n", worst);
  printf("large: worst relative residual %g\n", large_trials(rng, 300));
}
