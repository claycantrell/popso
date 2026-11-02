// Tests for the LU factorization: hyper-sparse solves against the dense loops on random sparse
// bases (with Forrest-Tomlin / product-form updates, dense kernels, singular bases, renumbered
// positions).

#include <cmath>
#include <random>
#include <vector>

#include "milp/lp/lu.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

using Col = std::vector<std::pair<Index, double>>;

// LP-like basis: unit columns, short structural columns, optionally a dense block (dense kernel
// path of the factorization), a duplicated and an empty column (replaced as singular).
std::vector<Col> random_basis(std::mt19937& rng, Index m, bool dense_block) {
  std::vector<Col> cols(m);
  const int slack_pct = static_cast<int>(rng() % 90), len = 1 + static_cast<int>(rng() % 6);
  for (Index c = 0; c < m; ++c) {
    if (static_cast<int>(rng() % 100) < slack_pct) {
      cols[c].push_back({static_cast<Index>(rng() % m), rng() % 3 ? 1.0 : -2.5});
      continue;
    }
    std::vector<char> used(m, 0);
    for (int k = 0; k < len; ++k) {
      const Index r = static_cast<Index>(rng() % m);
      if (used[r]) continue;
      used[r] = 1;
      cols[c].push_back({r, static_cast<double>(static_cast<int>(rng() % 19) - 9) + ((rng() % 4) ? 0.0 : 0.37)});
    }
  }
  if (dense_block && m >= 120) {
    const Index k = 60 + static_cast<Index>(rng() % 20), r0 = static_cast<Index>(rng() % (m - k));
    for (Index c = 0; c < k; ++c) {
      Col col;
      for (Index i = 0; i < k; ++i)
        if (rng() % 10 < 7 || i == c) col.push_back({r0 + i, 0.5 + static_cast<double>(rng() % 100) / 37.0});
      cols[m - 1 - c] = col;
    }
  }
  if (rng() % 4 == 0 && m > 4) cols[2] = cols[1];
  if (rng() % 5 == 0 && m > 4) cols[3].clear();
  return cols;
}

void factorize(LuFactor& lu, std::vector<Col>& cols) {
  const Index m = static_cast<Index>(cols.size());
  std::vector<Index> st{0}, ix;
  std::vector<double> vv;
  for (const Col& c : cols) {
    for (auto [r, v] : c) ix.push_back(r), vv.push_back(v);
    st.push_back(static_cast<Index>(ix.size()));
  }
  std::vector<std::pair<Index, Index>> rep;
  lu.factorize(m, st, ix, vv, &rep);
  for (auto [p, r] : rep) cols[p] = {{r, 1.0}};
}

HVector sparse_rhs(std::mt19937& rng, Index m, int nnz) {
  HVector h;
  h.resize(m);
  for (int k = 0; k < nnz; ++k) {
    const Index i = static_cast<Index>(rng() % m);
    if (h.dense[i] == 0.0) h.nz[h.count++] = i;
    h.dense[i] += static_cast<double>(static_cast<int>(rng() % 9) - 4) + 0.25;
    if (h.dense[i] == 0.0) h.dense[i] = 1.0;
  }
  return h;
}

// a and b agree to rounding and both index lists are exact (every nonzero listed once).
bool same(const HVector& a, const HVector& b, double* err) {
  const Index m = static_cast<Index>(a.dense.size());
  double mx = 0.0, e = 0.0;
  for (Index i = 0; i < m; ++i) {
    mx = std::max(mx, std::abs(a.dense[i]));
    e = std::max(e, std::abs(a.dense[i] - b.dense[i]));
  }
  *err = e / (1.0 + mx);
  for (const HVector* v : {&a, &b}) {
    std::vector<int> listed(m, 0);
    for (Index k = 0; k < v->count; ++k) ++listed[v->nz[k]];
    for (Index i = 0; i < m; ++i) {
      if (listed[i] > 1) return false;
      if (v->dense[i] != 0.0 && !listed[i]) return false;
    }
  }
  return *err <= 1e-11;
}

}  // namespace

TEST(lu_hyper_sparse_solves_match_dense) {
  std::mt19937 rng(17);
  double worst = 0.0;
  int compared = 0;
  for (int trial = 0; trial < 120; ++trial) {
    const Index m = 20 + static_cast<Index>(rng() % (trial % 3 == 0 ? 700 : 250));
    std::vector<Col> cols = random_basis(rng, m, trial % 4 == 1);
    LuFactor lu;
    if (trial % 5 == 4) lu.forrest_tomlin = false;
    factorize(lu, cols);
    if (trial % 2 == 0) {
      std::vector<Index> old_pos;
      lu.renumber_positions(old_pos);
      std::vector<Col> perm(m);
      for (Index p = 0; p < m; ++p) perm[p] = cols[old_pos[p]];
      cols.swap(perm);
    }
    LuFactor dense = lu, hyper = lu;
    dense.hyper_mode = 1;
    hyper.hyper_mode = 2;
    const int updates = static_cast<int>(rng() % 60);
    for (int u = 0; u <= updates; ++u) {
      // Compare solves on sparse and denser right-hand sides.
      for (int s = 0; s < 3; ++s) {
        const int nnz = s == 2 ? static_cast<int>(m / 3) : 1 + static_cast<int>(rng() % 4);
        HVector f1 = sparse_rhs(rng, m, nnz), f2 = f1, b1 = sparse_rhs(rng, m, nnz), b2 = b1;
        dense.ftran(f1);
        hyper.ftran(f2);
        dense.btran(b1);
        hyper.btran(b2);
        double e1, e2;
        CHECK(same(f1, f2, &e1));
        CHECK(same(b1, b2, &e2));
        worst = std::max({worst, e1, e2});
        compared += 2;
      }
      if (u == updates) break;
      // Replace a column: the ftran of the new column (hyper-sparse or dense, alternating) is the
      // last ftran before the update of both copies, as the update requires.
      Col a;
      for (int k = 0; k < 1 + static_cast<int>(rng() % 4); ++k)
        a.push_back({static_cast<Index>(rng() % m), static_cast<double>(static_cast<int>(rng() % 7) - 3)});
      HVector h1, h2;
      h1.resize(m);
      for (auto [r, v] : a) {
        if (h1.dense[r] == 0.0) h1.nz[h1.count++] = r;
        h1.dense[r] += v;
      }
      h2 = h1;
      dense.ftran(h1);
      hyper.ftran(h2);
      double e;
      CHECK(same(h1, h2, &e));
      const Index p = static_cast<Index>(rng() % m);
      double amax = 0.0;
      for (Index i = 0; i < m; ++i) amax = std::max(amax, std::abs(h1.dense[i]));
      if (std::abs(h1.dense[p]) < std::max(0.1 * amax, 1e-3)) {
        // Unsuitable pivot: still keep the copies in step (their next update needs a fresh ftran).
        continue;
      }
      dense.update(p, h1);
      hyper.update(p, h2);
      Col ac;
      for (Index r = 0; r < m; ++r) {
        double v = 0.0;
        for (auto [rr, vv] : a)
          if (rr == r) v += vv;
        if (v != 0.0) ac.push_back({r, v});
      }
      cols[p] = ac;
    }
    // Residual of the hyper-sparse solve against the basis: B x = b.
    HVector b = sparse_rhs(rng, m, 3), x = b;
    hyper.ftran(x);
    std::vector<double> bx(m, 0.0);
    for (Index c = 0; c < m; ++c)
      for (auto [r, v] : cols[c]) bx[r] += v * x.dense[c];
    double res = 0.0, scale = 1.0;
    for (Index i = 0; i < m; ++i) {
      res = std::max(res, std::abs(bx[i] - b.dense[i]));
      scale = std::max(scale, std::abs(b.dense[i]));
    }
    CHECK(res <= 1e-7 * scale);
  }
  CHECK(compared > 1000);
  CHECK(worst <= 1e-11);
}

TEST(lu_multi_ftran_mixed_modes) {
  // ftran(a, b, c) with some right-hand sides solved hyper-sparsely and the others in one dense pass.
  std::mt19937 rng(5);
  for (int trial = 0; trial < 40; ++trial) {
    const Index m = 50 + static_cast<Index>(rng() % 400);
    std::vector<Col> cols = random_basis(rng, m, trial % 3 == 0);
    LuFactor lu;
    factorize(lu, cols);
    HVector a = sparse_rhs(rng, m, 2), b = sparse_rhs(rng, m, static_cast<int>(m / 2)), c = sparse_rhs(rng, m, 3);
    HVector a1 = a, b1 = b, c1 = c;
    lu.ftran(a, b, c, 0.0, 0.9, 0.0);  // a, c hyper-sparse, b dense
    lu.hyper_mode = 1;
    lu.ftran(a1);
    lu.ftran(b1);
    lu.ftran(c1);
    double e;
    CHECK(same(a, a1, &e));
    CHECK(same(b, b1, &e));
    CHECK(same(c, c1, &e));
  }
}
