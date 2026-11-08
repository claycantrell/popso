#include "milp/conic/socp.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>

#include "milp/core/interrupt.hpp"
#include "milp/core/timer.hpp"
#include "milp/lp/cholesky.hpp"

namespace milp {

const char* to_string(ConicStatus s) {
  switch (s) {
    case ConicStatus::Optimal: return "optimal";
    case ConicStatus::PrimalInfeasible: return "primal infeasible";
    case ConicStatus::DualInfeasible: return "dual infeasible";
    case ConicStatus::IterationLimit: return "iteration limit";
    case ConicStatus::TimeLimit: return "time limit";
    case ConicStatus::NumericalError: return "numerical error";
  }
  return "?";
}

namespace {

double dot(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0.0;
  for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

double norm2(const std::vector<double>& a) { return std::sqrt(dot(a, a)); }

// out = M v (M CSC with n columns, out sized to its rows), or out += sign * M' v.
void mul(const SparseMatrix& M, const std::vector<double>& v, std::vector<double>& out) {
  out.assign(static_cast<size_t>(M.num_minor), 0.0);
  for (Index j = 0; j < M.num_major; ++j) {
    const double vj = v[j];
    if (vj == 0.0) continue;
    for (Index k = M.begin(j); k < M.end(j); ++k) out[M.index[k]] += M.value[k] * vj;
  }
}
// out += S v for S symmetric given by its lower triangle L (CSC, with the diagonal).
void sym_mul_add(const SparseMatrix& L, const std::vector<double>& v, std::vector<double>& out) {
  for (Index j = 0; j < L.num_major; ++j) {
    for (Index k = L.begin(j); k < L.end(j); ++k) {
      const Index i = L.index[k];
      out[i] += L.value[k] * v[j];
      if (i != j) out[j] += L.value[k] * v[i];
    }
  }
}

void mul_t_add(const SparseMatrix& M, const std::vector<double>& v, double sign, std::vector<double>& out) {
  for (Index j = 0; j < M.num_major; ++j) {
    double s = 0.0;
    for (Index k = M.begin(j); k < M.end(j); ++k) s += M.value[k] * v[M.index[k]];
    out[j] += sign * s;
  }
}

// The cone K = R+^l x Q^{q_1} x ... and the Nesterov-Todd scaling W (symmetric) at a pair (s, z) of its interior:
// on the orthant W = diag(sqrt(s / z)); on a cone W = eta (2 w w' - J)^{1/2}, i.e.
//   W = eta [ w0  w1' ; w1  I + w1 w1' / (1 + w0) ],  W^2 = eta^2 (2 w w' - J),  w0^2 - ||w1||^2 = 1,
// with W z = W^-1 s = lambda (Vandenberghe, "The CVXOPT linear and quadratic cone program solvers", 2010).
struct Cone {
  Index m = 0, l = 0;
  std::vector<Index> first, size;  // second-order cones: first row, rows
  std::vector<double> wl;          // orthant scaling sqrt(s / z)
  std::vector<double> eta;         // per cone
  std::vector<double> w;           // per cone row: the scaling point (w0, w1)
  Index degree() const { return l + static_cast<Index>(first.size()); }

  void init(Index l_, const std::vector<Index>& cones) {
    l = l_;
    first.clear(), size.clear();
    Index r = l;
    for (Index q : cones) first.push_back(r), size.push_back(q), r += q;
    m = r;
    wl.assign(static_cast<size_t>(l), 1.0);
    eta.assign(first.size(), 1.0);
    w.assign(static_cast<size_t>(m), 0.0);
    for (size_t k = 0; k < first.size(); ++k) w[first[k]] = 1.0;
  }

  // s0^2 - ||s1||^2 of cone k, as (s0 - ||s1||) (s0 + ||s1||) (less cancellation near the boundary)
  double det(const double* v, size_t k) const {
    const Index f = first[k], q = size[k];
    double t = 0.0;
    for (Index i = 1; i < q; ++i) t += v[f + i] * v[f + i];
    const double nt = std::sqrt(t);
    return (v[f] - nt) * (v[f] + nt);
  }

  // v + a d strictly inside K (each cone's s0 - ||s1|| > 0, the orthant's entries > 0).
  bool inside(const std::vector<double>& v, const std::vector<double>& d, double a) const {
    for (Index i = 0; i < l; ++i) {
      if (!(v[i] + a * d[i] > 0.0)) return false;
    }
    for (size_t k = 0; k < first.size(); ++k) {
      const Index f = first[k], q = size[k];
      double t = 0.0;
      for (Index i = 1; i < q; ++i) t += (v[f + i] + a * d[f + i]) * (v[f + i] + a * d[f + i]);
      if (!(v[f] + a * d[f] - std::sqrt(t) > 0.0)) return false;
    }
    return true;
  }

  // The unit e of K.
  void unit(std::vector<double>& e, double scale) const {
    e.assign(static_cast<size_t>(m), 0.0);
    for (Index i = 0; i < l; ++i) e[i] = scale;
    for (Index f : first) e[f] = scale;
  }

  // Moves v into the interior: v + (1 + a) e where a = the smallest shift that reaches the boundary, if v is not
  // already well inside (Vandenberghe's initialization).
  void shift_inside(std::vector<double>& v) const {
    double a = -kInf;
    for (Index i = 0; i < l; ++i) a = std::max(a, -v[i]);
    for (size_t k = 0; k < first.size(); ++k) {
      const Index f = first[k], q = size[k];
      double t = 0.0;
      for (Index i = 1; i < q; ++i) t += v[f + i] * v[f + i];
      a = std::max(a, std::sqrt(t) - v[f]);
    }
    if (a < 0.0 && a > -kInf) return;  // strictly inside
    if (a == -kInf) a = 0.0;
    const double add = 1.0 + std::max(0.0, a);
    for (Index i = 0; i < l; ++i) v[i] += add;
    for (Index f : first) v[f] += add;
  }

  // Largest step a >= 0 with v + a d in K (v in its interior); kInf if unbounded.
  double max_step(const std::vector<double>& v, const std::vector<double>& d) const {
    double a = kInf;
    for (Index i = 0; i < l; ++i) {
      if (d[i] < 0.0) a = std::min(a, -v[i] / d[i]);
    }
    for (size_t k = 0; k < first.size(); ++k) {
      const Index f = first[k], q = size[k];
      // f(t) = C + 2 B t + A t^2 >= 0 (and the cone's own half): the first positive root, if any.
      double A = d[f] * d[f], B = v[f] * d[f], C = v[f] * v[f];
      for (Index i = 1; i < q; ++i) A -= d[f + i] * d[f + i], B -= v[f + i] * d[f + i], C -= v[f + i] * v[f + i];
      if (C <= 0.0 || v[f] <= 0.0) return 0.0;
      const double disc = B * B - A * C;
      if (disc < 0.0 && A >= 0.0) continue;  // never reaches the boundary
      const double den = -B + std::sqrt(std::max(0.0, disc));
      if (den > 0.0) a = std::min(a, C / den);
    }
    return a;
  }

  // NT scaling at (s, z); false if they are not interior.
  bool scale(const std::vector<double>& s, const std::vector<double>& z) {
    for (Index i = 0; i < l; ++i) {
      if (!(s[i] > 0.0 && z[i] > 0.0)) return false;
      wl[i] = std::sqrt(s[i] / z[i]);
    }
    for (size_t k = 0; k < first.size(); ++k) {
      const Index f = first[k], q = size[k];
      const double ds = det(s.data(), k), dz = det(z.data(), k);
      if (!(ds > 0.0 && dz > 0.0 && s[f] > 0.0 && z[f] > 0.0)) return false;
      const double sr = std::sqrt(ds), zr = std::sqrt(dz);
      double sz = 0.0;  // s_bar' z_bar
      for (Index i = 0; i < q; ++i) sz += (s[f + i] / sr) * (z[f + i] / zr);
      const double gamma = std::sqrt(std::max(0.5 * (1.0 + sz), 1e-300));
      w[f] = (s[f] / sr + z[f] / zr) / (2.0 * gamma);
      for (Index i = 1; i < q; ++i) w[f + i] = (s[f + i] / sr - z[f + i] / zr) / (2.0 * gamma);
      // Renormalize w0 so that w0^2 - ||w1||^2 = 1 (round-off).
      double t = 0.0;
      for (Index i = 1; i < q; ++i) t += w[f + i] * w[f + i];
      w[f] = std::sqrt(1.0 + t);
      eta[k] = std::sqrt(sr / zr);
    }
    return true;
  }

  // out = W v, or W^-1 v (inverse).
  void apply(const std::vector<double>& v, std::vector<double>& out, bool inverse) const {
    out.resize(static_cast<size_t>(m));
    for (Index i = 0; i < l; ++i) out[i] = inverse ? v[i] / wl[i] : v[i] * wl[i];
    for (size_t k = 0; k < first.size(); ++k) {
      const Index f = first[k], q = size[k];
      const double w0 = w[f];
      double zeta = 0.0;
      for (Index i = 1; i < q; ++i) zeta += w[f + i] * v[f + i];
      const double e = inverse ? 1.0 / eta[k] : eta[k];
      const double v0 = v[f];
      out[f] = e * (w0 * v0 + (inverse ? -zeta : zeta));
      const double coef = (inverse ? -v0 : v0) + zeta / (1.0 + w0);
      for (Index i = 1; i < q; ++i) out[f + i] = e * (v[f + i] + coef * w[f + i]);
    }
  }

  // out = u o v (Jordan product: componentwise on the orthant, (u'v, u0 v1 + v0 u1) on a cone).
  void jordan(const std::vector<double>& u, const std::vector<double>& v, std::vector<double>& out) const {
    out.resize(static_cast<size_t>(m));
    for (Index i = 0; i < l; ++i) out[i] = u[i] * v[i];
    for (size_t k = 0; k < first.size(); ++k) {
      const Index f = first[k], q = size[k];
      double t = 0.0;
      for (Index i = 0; i < q; ++i) t += u[f + i] * v[f + i];
      for (Index i = 1; i < q; ++i) out[f + i] = u[f] * v[f + i] + v[f] * u[f + i];
      out[f] = t;
    }
  }

  // out = lambda \ d: the u with lambda o u = d.
  void jordan_inv(const std::vector<double>& lam, const std::vector<double>& d, std::vector<double>& out) const {
    out.resize(static_cast<size_t>(m));
    for (Index i = 0; i < l; ++i) out[i] = d[i] / lam[i];
    for (size_t k = 0; k < first.size(); ++k) {
      const Index f = first[k], q = size[k];
      const double rho = det(lam.data(), k);
      double ld = 0.0;
      for (Index i = 1; i < q; ++i) ld += lam[f + i] * d[f + i];
      const double u0 = (lam[f] * d[f] - ld) / rho;
      out[f] = u0;
      for (Index i = 1; i < q; ++i) out[f + i] = (d[f + i] - u0 * lam[f + i]) / lam[f];
    }
  }

  // out = W^2 v.
  void apply_w2(const std::vector<double>& v, std::vector<double>& out) const {
    out.resize(static_cast<size_t>(m));
    for (Index i = 0; i < l; ++i) out[i] = wl[i] * wl[i] * v[i];
    for (size_t k = 0; k < first.size(); ++k) {
      const Index f = first[k], q = size[k];
      double wv = 0.0;  // w' v
      for (Index i = 0; i < q; ++i) wv += w[f + i] * v[f + i];
      const double e2 = eta[k] * eta[k];
      out[f] = e2 * (2.0 * w[f] * wv - v[f]);
      for (Index i = 1; i < q; ++i) out[f + i] = e2 * (2.0 * w[f + i] * wv + v[f + i]);
    }
  }
};

// The KKT matrix (lower triangle, CSC) with its pattern fixed and the value positions to update. Columns: x, y, z,
// then two per large cone (a_k, b_k). A small cone (at most kExpandAbove rows) has its -W^2 as a dense block. A large
// one uses the sparse expansion of its scaling (as in ECOS),
//   W^2 = eta^2 (2 w w' - J) = eta^2 (D + u u' - v v'),  D = diag(d1, 1, ..., 1),  u = (u0, u1 w1),  v = (0, v1 w1),
// with the column a_k = -eta u (pivot +1) and b_k = eta v (pivot -1) against the cone's z rows (diagonal -eta^2 D):
// eliminating them gives back -W^2, at O(q) entries instead of q^2 / 2. With v1^2 = (2 / (2t + 1) + 1 / t) / 2
// (t = |w1|^2 = w0^2 - 1, > 1/4; else 3), u1 = sqrt(2 + v1^2), u0 = 2 w0 / u1 and d1 = 1 / (2 t (2 + v1^2)) (else
// (v1^2 (2t + 1) - 2) / (2 + v1^2)), d1 > 0 and D - v v' is positive definite (v1^2 t < 1): the matrix stays
// quasi-definite (x, a positive; y, z, b negative). But d1 ~ 1 / (4 w0^2) is tiny near the cone's boundary, so the
// expansion is numerically weaker than the dense block (QPLIB_2482: refinement failed at iteration 15 with it, 4-row
// cones): only for large cones.
struct Kkt {
  static constexpr Index kExpandAbove = 32;
  Index n = 0, p = 0, m = 0, N = 0;  // N: the factored system (x, y, the kept z rows, the expansions' columns)
  std::vector<Index> start, index;
  std::vector<double> value;
  std::vector<size_t> diag;      // per factored column: position of its diagonal
  std::vector<Index> zpos;       // per z row: its factored column, -1 for an eliminated singleton row
  std::vector<size_t> cone_pos;  // dense cone: position of its block's first entry; expanded: unused
  std::vector<Index> ab;         // expanded cone: its column a_k (b_k = a_k + 1); dense: -1
  std::vector<size_t> cone_off;  // per z row of an expanded cone: position of its a-entry (the b-entry follows)
  // Orthant rows with one entry g in G (bounds, mostly) are eliminated: dz = (g dx_j - r_z) / (W^2 + delta), which adds
  // g^2 / (W^2 + delta) to x_j's pivot and g r_z / (W^2 + delta) to its right side (QPLIB_8602: 69104 bound rows made
  // the KKT 156639 wide, a 3.5e7-nonzero factor).
  std::vector<Index> sg_col;     // per z row: the singleton's column (-1: not eliminated)
  std::vector<double> sg_val, sg_w2;
  std::vector<int8_t> signs;
  std::vector<double> pdiag;  // P's diagonal (added to the x pivots)
  SparseCholesky chol;
  double delta = 1e-8;  // static regularization

  bool build(const SparseMatrix& P, const SparseMatrix& A, const SparseMatrix& G, const Cone& K, int threads) {
    n = A.num_major, p = A.num_minor, m = G.num_minor;
    const Index nc = static_cast<Index>(K.first.size());
    sg_col.assign(static_cast<size_t>(m), -1);
    sg_val.assign(static_cast<size_t>(m), 0.0);
    sg_w2.assign(static_cast<size_t>(m), 1.0);
    {
      std::vector<Index> cnt(static_cast<size_t>(m), 0);
      for (Index j = 0; j < n; ++j) {
        for (Index k = G.begin(j); k < G.end(j); ++k) {
          const Index i = G.index[k];
          if (i < K.l && G.value[k] != 0.0) ++cnt[i], sg_col[i] = j, sg_val[i] = G.value[k];
        }
      }
      for (Index i = 0; i < m; ++i) {
        if (i >= K.l || cnt[i] != 1) sg_col[i] = -1;
      }
    }
    zpos.assign(static_cast<size_t>(m), -1);
    N = n + p;
    for (Index i = 0; i < m; ++i) {
      if (sg_col[i] < 0) zpos[i] = N++;
    }
    ab.assign(static_cast<size_t>(nc), -1);
    for (Index k = 0; k < nc; ++k) {
      if (K.size[k] > kExpandAbove) ab[k] = N, N += 2;
    }
    start.assign(1, 0), index.clear();
    diag.assign(static_cast<size_t>(N), 0);
    cone_pos.assign(static_cast<size_t>(nc), 0);
    cone_off.assign(static_cast<size_t>(m), 0);
    std::vector<Index> cone_of(static_cast<size_t>(m), -1);
    for (Index k = 0; k < nc; ++k) {
      for (Index i = 0; i < K.size[k]; ++i) cone_of[K.first[k] + i] = k;
    }
    const bool hasP = P.num_major == n;
    pdiag.assign(static_cast<size_t>(n), 0.0);
    for (Index j = 0; j < n; ++j) {
      diag[j] = index.size();
      index.push_back(j);
      if (hasP) {
        for (Index k = P.begin(j); k < P.end(j); ++k) {
          if (P.index[k] > j) index.push_back(P.index[k]);
          else if (P.index[k] == j) pdiag[j] += P.value[k];
        }
      }
      for (Index k = A.begin(j); k < A.end(j); ++k) index.push_back(n + A.index[k]);
      for (Index k = G.begin(j); k < G.end(j); ++k) {
        if (zpos[G.index[k]] >= 0) index.push_back(zpos[G.index[k]]);
      }
      start.push_back(static_cast<Index>(index.size()));
    }
    for (Index i = 0; i < p; ++i) {
      diag[n + i] = index.size();
      index.push_back(n + i);
      start.push_back(static_cast<Index>(index.size()));
    }
    for (Index r = 0; r < m; ++r) {
      if (zpos[r] < 0) continue;
      diag[zpos[r]] = index.size();
      index.push_back(zpos[r]);
      if (const Index k = cone_of[r]; k >= 0) {
        if (ab[k] >= 0) {
          cone_off[r] = index.size();
          index.push_back(ab[k]);
          index.push_back(ab[k] + 1);
        } else {
          if (r == K.first[k]) cone_pos[k] = diag[zpos[r]];
          for (Index r2 = r + 1; r2 < K.first[k] + K.size[k]; ++r2) index.push_back(zpos[r2]);
        }
      }
      start.push_back(static_cast<Index>(index.size()));
    }
    for (Index k = 0; k < nc; ++k) {
      for (Index c = ab[k]; ab[k] >= 0 && c < ab[k] + 2; ++c) {
        diag[c] = index.size();
        index.push_back(c);
        start.push_back(static_cast<Index>(index.size()));
      }
    }
    value.assign(index.size(), 0.0);
    // The P, A and G values never change.
    for (Index j = 0; j < n; ++j) {
      size_t pos = diag[j] + 1;
      if (hasP) {
        for (Index k = P.begin(j); k < P.end(j); ++k) {
          if (P.index[k] > j) value[pos++] = P.value[k];
        }
      }
      for (Index k = A.begin(j); k < A.end(j); ++k) value[pos++] = A.value[k];
      for (Index k = G.begin(j); k < G.end(j); ++k) {
        if (zpos[G.index[k]] >= 0) value[pos++] = G.value[k];
      }
    }
    signs.assign(static_cast<size_t>(N), -1);
    for (Index j = 0; j < n; ++j) signs[j] = 1;
    for (Index k = 0; k < nc; ++k) {
      if (ab[k] >= 0) signs[ab[k]] = 1;
    }
    chol.threads = threads;
    chol.interruptible = true;
    chol.dynamic_reg = 1e-8;
    chol.pivot_tol = 1e-13;
    return chol.analyze(N, start, index, nullptr, &signs);
  }

  // The scaling (W = I when ident) and the regularization; then the factorization.
  bool factor(const Cone& K, bool ident) {
    for (Index j = 0; j < n; ++j) value[diag[j]] = pdiag[j] + delta;
    for (Index i = 0; i < p; ++i) value[diag[n + i]] = -delta;
    for (Index i = 0; i < K.l; ++i) {
      const double w2 = (ident ? 1.0 : K.wl[i] * K.wl[i]) + delta;
      if (sg_col[i] >= 0) {
        sg_w2[i] = w2;
        value[diag[sg_col[i]]] += sg_val[i] * sg_val[i] / w2;
      } else {
        value[diag[zpos[i]]] = -w2;
      }
    }
    for (size_t k = 0; k < K.first.size(); ++k) {
      const Index f = K.first[k], q = K.size[k];
      const double eta = ident ? 1.0 : K.eta[k], e2 = eta * eta;
      if (ab[k] < 0) {
        // Dense block: W^2 = eta^2 (2 w w' - J), column-major lower part.
        size_t pos = cone_pos[k];
        for (Index a = 0; a < q; ++a) {
          for (Index b = a; b < q; ++b) {
            double v = ident ? (a == b ? 1.0 : 0.0) : 2.0 * K.w[f + a] * K.w[f + b];
            if (!ident && a == b) v += (a == 0 ? -1.0 : 1.0);
            value[pos++] = -e2 * v - (a == b ? delta : 0.0);
          }
        }
        continue;
      }
      value[diag[ab[k]]] = 1.0;
      value[diag[ab[k] + 1]] = -1.0;
      if (ident) {
        for (Index i = 0; i < q; ++i) {
          value[diag[zpos[f + i]]] = -1.0 - delta;
          value[cone_off[f + i]] = value[cone_off[f + i] + 1] = 0.0;
        }
        continue;
      }
      const double w0 = K.w[f];
      double t = 0.0;
      for (Index i = 1; i < q; ++i) t += K.w[f + i] * K.w[f + i];
      const double v1sq = t > 0.25 ? 0.5 * (2.0 / (2.0 * t + 1.0) + 1.0 / t) : 3.0;
      const double u1 = std::sqrt(2.0 + v1sq), u0 = 2.0 * w0 / u1, v1 = std::sqrt(v1sq);
      const double d1 = t > 0.25 ? 1.0 / (2.0 * t * (2.0 + v1sq)) : (v1sq * (2.0 * t + 1.0) - 2.0) / (2.0 + v1sq);
      value[diag[zpos[f]]] = -e2 * d1 - delta;
      value[cone_off[f]] = -eta * u0;
      value[cone_off[f] + 1] = 0.0;
      for (Index i = 1; i < q; ++i) {
        value[diag[zpos[f + i]]] = -e2 - delta;
        value[cone_off[f + i]] = -eta * u1 * K.w[f + i];
        value[cone_off[f + i] + 1] = eta * v1 * K.w[f + i];
      }
    }
    return chol.factorize(value);
  }

  // Solves the regularized system in place for v = (x, y, z) (n + p + m): the singleton rows folded into x, the rest
  // factored (the expansions' columns with zero right sides), then the singletons' dz recovered.
  void solve(std::vector<double>& v, std::vector<double>& work) const {
    work.assign(static_cast<size_t>(N), 0.0);
    for (Index j = 0; j < n + p; ++j) work[j] = v[j];
    for (Index i = 0; i < m; ++i) {
      if (zpos[i] >= 0) work[zpos[i]] = v[n + p + i];
      else work[sg_col[i]] += sg_val[i] * v[n + p + i] / sg_w2[i];
    }
    chol.solve(work);
    for (Index i = 0; i < m; ++i) {
      v[n + p + i] = zpos[i] >= 0 ? work[zpos[i]] : (sg_val[i] * work[sg_col[i]] - v[n + p + i]) / sg_w2[i];
    }
    for (Index j = 0; j < n + p; ++j) v[j] = work[j];
  }
};

}  // namespace

ConicResult solve_conic(const ConicProblem& prob, const ConicOptions& opt) {
  Timer timer;
  ConicResult res;
  const Index n = prob.n, p = prob.A.num_minor, m = prob.G.num_minor;
  const Logger quiet(0);
  const Logger& log = opt.log ? *opt.log : quiet;
  Cone K;
  K.init(prob.l, prob.cones);
  if (K.m != m) {
    log.log(1, "conic: the cones cover %d of %d rows\n", K.m, m);
    return res;
  }

  // ---- Equilibration (Ruiz): A~ = E A D, G~ = F G D, c~ = D c, b~ = E b, h~ = F h.
  SparseMatrix A = prob.A, G = prob.G, P = prob.P;
  const bool hasP = P.num_major == n && P.nnz() > 0;
  if (!hasP) P = SparseMatrix{};
  std::vector<double> D(static_cast<size_t>(n), 1.0), E(static_cast<size_t>(p), 1.0), F(static_cast<size_t>(m), 1.0);
  for (int pass = 0; pass < opt.equilibration; ++pass) {
    std::vector<double> cmax(static_cast<size_t>(n), 0.0), amax(static_cast<size_t>(p), 0.0), gmax(static_cast<size_t>(m), 0.0);
    for (Index j = 0; j < n; ++j) {
      for (Index k = A.begin(j); k < A.end(j); ++k) {
        const double a = std::abs(A.value[k]);
        cmax[j] = std::max(cmax[j], a), amax[A.index[k]] = std::max(amax[A.index[k]], a);
      }
      for (Index k = G.begin(j); k < G.end(j); ++k) {
        const double a = std::abs(G.value[k]);
        cmax[j] = std::max(cmax[j], a), gmax[G.index[k]] = std::max(gmax[G.index[k]], a);
      }
      if (hasP) {
        for (Index k = P.begin(j); k < P.end(j); ++k) {
          const double a = std::abs(P.value[k]);
          cmax[j] = std::max(cmax[j], a), cmax[P.index[k]] = std::max(cmax[P.index[k]], a);
        }
      }
    }
    for (size_t k = 0; k < K.first.size(); ++k) {  // one factor per cone
      double t = 0.0;
      for (Index i = 0; i < K.size[k]; ++i) t = std::max(t, gmax[K.first[k] + i]);
      for (Index i = 0; i < K.size[k]; ++i) gmax[K.first[k] + i] = t;
    }
    auto f = [](double v) { return v > 0.0 ? 1.0 / std::sqrt(v) : 1.0; };
    for (Index j = 0; j < n; ++j) {
      const double dj = f(cmax[j]);
      D[j] *= dj;
      for (Index k = A.begin(j); k < A.end(j); ++k) A.value[k] *= dj;
      for (Index k = G.begin(j); k < G.end(j); ++k) G.value[k] *= dj;
    }
    if (hasP) {
      for (Index j = 0; j < n; ++j) {
        for (Index k = P.begin(j); k < P.end(j); ++k) P.value[k] *= f(cmax[j]) * f(cmax[P.index[k]]);
      }
    }
    for (Index i = 0; i < p; ++i) E[i] *= f(amax[i]);
    for (Index i = 0; i < m; ++i) F[i] *= f(gmax[i]);
    for (Index j = 0; j < n; ++j) {
      for (Index k = A.begin(j); k < A.end(j); ++k) A.value[k] *= f(amax[A.index[k]]);
      for (Index k = G.begin(j); k < G.end(j); ++k) G.value[k] *= f(gmax[G.index[k]]);
    }
  }
  std::vector<double> c(static_cast<size_t>(n)), b(static_cast<size_t>(p)), h(static_cast<size_t>(m));
  for (Index j = 0; j < n; ++j) c[j] = prob.c[j] * D[j];
  for (Index i = 0; i < p; ++i) b[i] = prob.b[i] * E[i];
  for (Index i = 0; i < m; ++i) h[i] = prob.h[i] * F[i];
  const double nb = std::max(1.0, norm2(b)), nh = std::max(1.0, norm2(h)), nc = std::max(1.0, norm2(c));

  Kkt kkt;
  if (!kkt.build(P, A, G, K, opt.threads)) {
    log.log(1, "conic: KKT analysis failed\n");
    return res;
  }
  log.log(1, "Conic barrier: %d columns, %d equality and %d cone rows (%d orthant, %zu cones), KKT %d x %d, factor %lld "
          "nonzeros, %.2e flops (%d bound rows eliminated)\n", n, p, m, K.l, K.first.size(), kkt.N, kkt.N,
          kkt.chol.factor_nnz(), kkt.chol.flops(), static_cast<int>(n + p + m + 2 * static_cast<Index>(std::count_if(
              kkt.ab.begin(), kkt.ab.end(), [](Index a) { return a >= 0; })) - kkt.N));
  const Index N = n + p + m;
  // K0 u (the unregularized KKT matrix with the current scaling; ident: W = I).
  auto kkt_mul = [&](const std::vector<double>& u, std::vector<double>& out, bool ident) {
    out.assign(static_cast<size_t>(N), 0.0);
    std::vector<double> ux(u.begin(), u.begin() + n), uy(u.begin() + n, u.begin() + n + p), uz(u.begin() + n + p, u.end());
    std::vector<double> t(static_cast<size_t>(n), 0.0), ay, gz, w2;
    mul_t_add(A, uy, 1.0, t);
    mul_t_add(G, uz, 1.0, t);
    if (hasP) sym_mul_add(P, ux, t);
    mul(A, ux, ay);
    mul(G, ux, gz);
    if (ident) {
      w2 = uz;
    } else {
      K.apply_w2(uz, w2);
    }
    for (Index j = 0; j < n; ++j) out[j] = t[j];
    for (Index i = 0; i < p; ++i) out[n + i] = ay[i];
    for (Index i = 0; i < m; ++i) out[n + p + i] = gz[i] - w2[i];
  };
  // Solves K0 u = r with the regularized factor and iterative refinement.
  // The iteration's worst residuals (relative): against the unregularized matrix, and against the regularized one (the
  // factor's own accuracy; the unregularized system can be singular: a free column in no row and not in P).
  double worst_solve = 0.0, worst_reg = 0.0;
  bool exact_guard = true;  // off once a singular system had to be accepted
  static const double kRegTol = [] {
    const char* e = std::getenv("MILP_CONIC_REGTOL");
    return e ? std::atof(e) : 1e-5;
  }();
  // Returns the final residual's max norm relative to 1 + |r|. Refinement keeps the best solution and stops when the
  // residual no longer shrinks (with regularized pivots it can diverge).
  // The factored system differs from (x, y, z): singleton rows eliminated, expansion columns added (Kkt::solve).
  std::vector<double> ext;
  auto ext_solve = [&](std::vector<double>& v) { kkt.solve(v, ext); };
  auto kkt_solve = [&](const std::vector<double>& r, std::vector<double>& u, bool ident) {
    u = r;
    ext_solve(u);
    std::vector<double> ku, d, best;
    const double nr = std::max(1.0, norm2(r));
    double e_best = kInf, ereg = kInf, ereg_best = kInf;
    int it = 0;
    for (;; ++it) {
      kkt_mul(u, ku, ident);
      d.resize(static_cast<size_t>(N));
      double e = 0.0;
      ereg = 0.0;
      for (Index i = 0; i < N; ++i) {
        d[i] = r[i] - ku[i], e = std::max(e, std::abs(d[i]));
        // The direct solve against the regularized matrix: the factor's own accuracy (K0 may be singular; after
        // refinement r - K_delta u ~ delta u says nothing).
        if (it == 0) ereg = std::max(ereg, std::abs(d[i] - (i < n ? kkt.delta : -kkt.delta) * u[i]));
      }
      if (it == 0) ereg_best = ereg;
      if (!(e < 0.5 * e_best)) {  // no progress (or not finite): keep the best
        if (!best.empty() && !(e < e_best)) u.swap(best);
        break;
      }
      e_best = e;
      if (e <= 1e-14 * nr || it >= 10) break;
      best = u;
      ext_solve(d);
      for (Index i = 0; i < N; ++i) u[i] += d[i];
    }
    log.log(3, "conic   refine %d steps, residual %.2e (regularized %.2e, rhs %.2e), regularized pivots %d\n", it, e_best,
            ereg_best, nr, static_cast<int>(kkt.chol.num_dropped()));
    if (!std::isfinite(e_best)) ereg_best = kInf;
    worst_solve = std::max(worst_solve, e_best / nr);
    worst_reg = std::max(worst_reg, ereg_best / nr);
    return e_best / nr;
  };

  // ---- Initialization (W = I): x, s from [0 A' G'; A 0 0; G 0 -I] (x, y, z) = (0, b, h), s = -z shifted into K;
  // y, z from the right side (-c, 0, 0), z shifted into K; tau = kappa = 1.
  if (!kkt.factor(K, true)) {
    log.log(1, "conic: initial factorization failed\n");
    return res;
  }
  std::vector<double> x(static_cast<size_t>(n)), y(static_cast<size_t>(p)), z(static_cast<size_t>(m)), s(static_cast<size_t>(m));
  {
    std::vector<double> r(static_cast<size_t>(N), 0.0), u;
    for (Index i = 0; i < p; ++i) r[n + i] = b[i];
    for (Index i = 0; i < m; ++i) r[n + p + i] = h[i];
    kkt_solve(r, u, true);
    for (Index j = 0; j < n; ++j) x[j] = u[j];
    for (Index i = 0; i < m; ++i) s[i] = -u[n + p + i];
    K.shift_inside(s);
    std::fill(r.begin(), r.end(), 0.0);
    for (Index j = 0; j < n; ++j) r[j] = -c[j];
    kkt_solve(r, u, true);
    for (Index i = 0; i < p; ++i) y[i] = u[n + i];
    for (Index i = 0; i < m; ++i) z[i] = u[n + p + i];
    K.shift_inside(z);
  }
  double tau = 1.0, kappa = 1.0;
  const double deg = static_cast<double>(K.degree());

  std::vector<double> rx, ry, rz, px, lam, tmp, tmp2, u1, u2, rhs(static_cast<size_t>(N));
  std::vector<double> dx(static_cast<size_t>(n)), dy(static_cast<size_t>(p)), dz(static_cast<size_t>(m)),
      ds(static_cast<size_t>(m));
  ConicStatus status = ConicStatus::IterationLimit;
  int iter = 0;
  double pres = kInf, dres = kInf, gap = kInf, pcost = 0.0, dcost = 0.0;
  // The best iterate within the reduced tolerances (1e-6, ECOS's "inaccurate" exit with tighter values): returned when
  // the method stops short of the full ones.
  struct Saved {
    std::vector<double> x, y, z, s;
    double tau = 0.0, kappa = 0.0, merit = kInf, pres = 0.0, dres = 0.0, gap = 0.0;
  } saved;
  const double kReduced = 1e-6;
  double delta_now = 1e-8;
  double stall_best = kInf;
  int stall_iter = 0;
  static const int kStallWindow = [] {  // MILP_CONIC_STALL: the stall window in iterations (default 10)
    const char* e = std::getenv("MILP_CONIC_STALL");
    return e ? std::atoi(e) : 10;
  }();
  for (;; ++iter) {
    // Residuals: rx = P x + A'y + G'z + c tau, ry = A x - b tau, rz = s + G x - h tau,
    // rtau = kappa + c'x + b'y + h'z + x'Px / tau.
    rx.assign(static_cast<size_t>(n), 0.0);
    mul_t_add(A, y, 1.0, rx);
    mul_t_add(G, z, 1.0, rx);
    px.assign(static_cast<size_t>(n), 0.0);
    if (hasP) sym_mul_add(P, x, px);
    const double xpx = dot(x, px);
    for (Index j = 0; j < n; ++j) rx[j] += px[j] + c[j] * tau;
    mul(A, x, ry);
    for (Index i = 0; i < p; ++i) ry[i] -= b[i] * tau;
    mul(G, x, rz);
    for (Index i = 0; i < m; ++i) rz[i] += s[i] - h[i] * tau;
    const double cx = dot(c, x), by = dot(b, y), hz = dot(h, z);
    const double rtau = kappa + cx + by + hz + xpx / tau;
    const double sz = dot(s, z);
    const double mu = (sz + tau * kappa) / (deg + 1.0);
    pres = std::max(norm2(ry) / nb, norm2(rz) / nh) / tau;
    dres = norm2(rx) / nc / tau;
    pcost = (0.5 * xpx / tau + cx) / tau, dcost = (-0.5 * xpx / tau - by - hz) / tau;
    gap = sz / (tau * tau);
    const double relgap = gap / std::max(1e-12, std::min(std::abs(pcost), std::abs(dcost)));
    // The objectives' difference too, relative to the objective with its constant: a constant that cancels most of
    // c'x + 1/2 x'Px makes the gap relative to the latter far smaller than the objective's own error (Maros-Meszaros
    // GOULDQP3: offset ~29648, s'z 2.8e-4, 1e-8 relative to |c'x| but 1.4e-4 of the objective 2.0628: "optimal" at
    // 2.062811 vs the optimum 2.062784; HS268: 2.6e-5, its optimum ~0).
    const double ogap = std::abs(pcost - dcost);
    const double orel =
        ogap / std::max(1.0, std::min(std::abs(pcost + prob.offset), std::abs(dcost + prob.offset)));
    log.log(3, "conic   ry %.2e rz %.2e tau %.2e\n", norm2(ry) / nb, norm2(rz) / nh, tau);
    log.log(2, "conic %3d  pcost %+.9e  dcost %+.9e  gap %.2e  pres %.2e  dres %.2e  k/t %.2e\n", iter,
            pcost + prob.offset, dcost + prob.offset, gap, pres, dres, kappa / tau);
    if (pres <= opt.feastol && dres <= opt.feastol && (gap <= opt.abstol || relgap <= opt.reltol) &&
        (ogap <= opt.abstol || orel <= opt.reltol)) {
      status = ConicStatus::Optimal;
      break;
    }
    if (pres <= kReduced && dres <= kReduced && (gap <= kReduced || relgap <= kReduced) &&
        (ogap <= kReduced || orel <= kReduced)) {
      const double merit = std::max({pres, dres, std::min(gap, relgap)});
      if (merit < saved.merit) saved = Saved{x, y, z, s, tau, kappa, merit, pres, dres, gap};
    }
    // Certificates (the embedding's tau -> 0): primal infeasibility (b'y + h'z < 0, A'y + G'z ~ 0), dual
    // infeasibility (c'x < 0, A x ~ 0, G x + s ~ 0, P x ~ 0).
    if (by + hz < 0.0) {
      std::vector<double> t(static_cast<size_t>(n), 0.0);
      mul_t_add(A, y, 1.0, t);
      mul_t_add(G, z, 1.0, t);
      if (norm2(t) / nc <= opt.feastol * -(by + hz) && tau < 1e-6 * kappa) {
        status = ConicStatus::PrimalInfeasible;
        break;
      }
    }
    if (cx < 0.0) {
      std::vector<double> ax, gx;
      mul(A, x, ax);
      mul(G, x, gx);
      for (Index i = 0; i < m; ++i) gx[i] += s[i];
      if (std::max({norm2(ax) / nb, norm2(gx) / nh, norm2(px) / nc}) <= opt.feastol * -cx && tau < 1e-6 * kappa) {
        status = ConicStatus::DualInfeasible;
        break;
      }
    }
    if (iter >= opt.max_iterations) {
      status = ConicStatus::IterationLimit;
      break;
    }
    if (timer.elapsed() > opt.time_limit || interrupted_or_expired()) {
      status = ConicStatus::TimeLimit;
      break;
    }
    // Stall: the worst of the residuals and gaps not halved in 10 iterations (QPLIB_10038: s'z -> 0 while the dual
    // residual crept from 1e-5 to 3e-6 over 30 iterations of 0.9 s; the QP barrier solves it in 5 s). The best saved
    // iterate if it is within the reduced tolerances, else a numerical error (the callers' fallbacks). Only on the way
    // to an optimum (kappa < tau): towards a certificate the residuals need not shrink.
    if (kappa < tau) {
      const double worst = std::max({pres, dres, std::min(gap, relgap), std::min(ogap, orel)});
      if (worst < 0.5 * stall_best) {
        stall_best = worst;
        stall_iter = iter;
      } else if (iter - stall_iter >= kStallWindow) {
        log.log(3, "conic   stop: stalled at %.2e since iteration %d\n", stall_best, stall_iter);
        status = ConicStatus::NumericalError;
        break;
      }
    }

    // ---- Scaling and factorization (static regularization delta_now: raised after a failed solve).
    kkt.delta = delta_now;
    worst_solve = 0.0, worst_reg = 0.0;
    const bool scaled = K.scale(s, z);
    // A factorization that overflows (large regularized pivots in a cone block): again with more static
    // regularization (refinement corrects for it), at most three times.
    bool factored = scaled && kkt.factor(K, false);
    const double delta0 = kkt.delta;
    for (int t = 0; scaled && !factored && t < 3; ++t) {
      kkt.delta *= 100.0;
      factored = kkt.factor(K, false);
      log.log(3, "conic   factorization again with delta %.1e: %s\n", kkt.delta, factored ? "ok" : "failed");
    }
    if (!factored) kkt.delta = delta0;
    if (!scaled || !factored) {
      double wmax = 0.0, emax = 0.0;
      for (double v : K.w) wmax = std::max(wmax, std::abs(v));
      for (double v : K.eta) emax = std::max(emax, v);
      log.log(3, "conic   stop: %s failed (max w %.2e, max eta %.2e)\n", scaled ? "factorization" : "scaling", wmax, emax);
      status = ConicStatus::NumericalError;
      break;
    }
    K.apply(z, lam, false);  // lambda = W z

    // The tau direction: K0 (x1, y1, z1) = (-c, b, h).
    for (Index j = 0; j < n; ++j) rhs[j] = -c[j];
    for (Index i = 0; i < p; ++i) rhs[n + i] = b[i];
    for (Index i = 0; i < m; ++i) rhs[n + p + i] = h[i];
    kkt_solve(rhs, u1, false);
    // With P: the tau row's x'Px / tau linearizes to (2 P x / tau)' dx - (x'Px / tau^2) dtau.
    double cx1 = 0.0, by1 = 0.0, hz1 = 0.0;
    for (Index j = 0; j < n; ++j) cx1 += (c[j] + 2.0 * px[j] / tau) * u1[j];
    for (Index i = 0; i < p; ++i) by1 += b[i] * u1[n + i];
    for (Index i = 0; i < m; ++i) hz1 += h[i] * u1[n + p + i];
    const double den = kappa / tau + xpx / (tau * tau) - cx1 - by1 - hz1;

    // One direction for the right sides d_x = -sx rx, d_y = sx ry, d_z = sx rz, d_tau = sx rtau, d_s, d_kappa:
    //   u = lambda \ d_s,  K0 (x2, y2, z2) = (d_x, -d_y, -d_z - W u),  Delta = (x2, y2, z2) + dtau (x1, y1, z1) with
    //   dtau = (d_tau + d_kappa / tau + c'x2 + b'y2 + h'z2) / den,  ds = W (u - W dz),  dkappa = (d_kappa - kappa dtau) / tau.
    double dtau = 0.0, dkappa = 0.0;
    std::vector<double> u, wu, wdz;
    auto direction = [&](double sx, const std::vector<double>& d_s, double d_kappa) {
      K.jordan_inv(lam, d_s, u);
      K.apply(u, wu, false);
      for (Index j = 0; j < n; ++j) rhs[j] = -sx * rx[j];
      for (Index i = 0; i < p; ++i) rhs[n + i] = -sx * ry[i];
      for (Index i = 0; i < m; ++i) rhs[n + p + i] = -sx * rz[i] - wu[i];
      kkt_solve(rhs, u2, false);
      double c2 = 0.0;
      for (Index j = 0; j < n; ++j) c2 += (c[j] + 2.0 * px[j] / tau) * u2[j];
      for (Index i = 0; i < p; ++i) c2 += b[i] * u2[n + i];
      for (Index i = 0; i < m; ++i) c2 += h[i] * u2[n + p + i];
      dtau = (sx * rtau + d_kappa / tau + c2) / den;
      for (Index j = 0; j < n; ++j) dx[j] = u2[j] + dtau * u1[j];
      for (Index i = 0; i < p; ++i) dy[i] = u2[n + i] + dtau * u1[n + i];
      for (Index i = 0; i < m; ++i) dz[i] = u2[n + p + i] + dtau * u1[n + p + i];
      K.apply(dz, wdz, false);
      tmp2.resize(static_cast<size_t>(m));
      for (Index i = 0; i < m; ++i) tmp2[i] = u[i] - wdz[i];
      K.apply(tmp2, ds, false);
      dkappa = (d_kappa - kappa * dtau) / tau;
    };
    auto step_to_boundary = [&]() {
      double a = std::min(K.max_step(s, ds), K.max_step(z, dz));
      if (dtau < 0.0) a = std::min(a, -tau / dtau);
      if (dkappa < 0.0) a = std::min(a, -kappa / dkappa);
      return a;
    };

    // Affine (predictor) direction: d_s = -lambda o lambda, d_kappa = -kappa tau.
    K.jordan(lam, lam, tmp);
    for (double& v : tmp) v = -v;
    direction(1.0, tmp, -kappa * tau);
    const double a_aff = std::min(1.0, step_to_boundary());
    const double sigma = std::clamp(std::pow(1.0 - a_aff, 3.0), 1e-4, 1.0);

    // Combined direction: d_s = -lambda o lambda - (W^-1 ds_a) o (W dz_a) + sigma mu e, d_kappa = -kappa tau -
    // dkappa_a dtau_a + sigma mu.
    {
      std::vector<double> wids, wdza, corr, e;
      K.apply(ds, wids, true);
      K.apply(dz, wdza, false);
      K.jordan(wids, wdza, corr);
      K.unit(e, sigma * mu);
      K.jordan(lam, lam, tmp);
      for (Index i = 0; i < m; ++i) tmp[i] = -tmp[i] - corr[i] + e[i];
      const double dk = -kappa * tau - dkappa * dtau + sigma * mu;
      direction(1.0 - sigma, tmp, dk);
    }
    if (!(worst_reg <= kRegTol) || (exact_guard && !(worst_solve <= 1e-3))) {
      // The Newton systems were not solved: again with more regularization, up to 1e-4. Then, if the factor solves its
      // own matrix, the unregularized one is singular (its least-squares-like solution is what the embedding needs
      // to detect unboundedness): go on and stop testing that; else stop (the best iterate, if within the reduced
      // tolerances).
      if (delta_now < 1e-4) {
        delta_now *= 100.0;
        log.log(3, "conic   solves failed (residual %.1e, regularized %.1e): again with delta %.0e\n", worst_solve,
                worst_reg, delta_now);
        continue;
      }
      if (!(worst_reg <= kRegTol)) {
        log.log(3, "conic   stop: solves failed (residual %.1e, regularized %.1e)\n", worst_solve, worst_reg);
        status = ConicStatus::NumericalError;
        break;
      }
      exact_guard = false;
    }
    delta_now = std::max(1e-8, delta_now / 10.0);
    double alpha = std::min(1.0, 0.99 * step_to_boundary());
    // The step's end must be strictly inside (the boundary step has round-off near it): back off if not.
    for (int t = 0; t < 30 && alpha > 0.0 && !(K.inside(s, ds, alpha) && K.inside(z, dz, alpha)); ++t) alpha *= 0.7;
    if (!(alpha > 1e-12) || !std::isfinite(alpha) || !(K.inside(s, ds, alpha) && K.inside(z, dz, alpha))) {
      log.log(3, "conic   stop: step %g (affine %g)\n", alpha, a_aff);
      status = ConicStatus::NumericalError;
      break;
    }
    for (Index j = 0; j < n; ++j) x[j] += alpha * dx[j];
    for (Index i = 0; i < p; ++i) y[i] += alpha * dy[i];
    for (Index i = 0; i < m; ++i) z[i] += alpha * dz[i], s[i] += alpha * ds[i];
    tau += alpha * dtau;
    kappa += alpha * dkappa;
    if (!(tau > 0.0) || !(kappa > 0.0)) {
      log.log(3, "conic   stop: tau %g kappa %g\n", tau, kappa);
      status = ConicStatus::NumericalError;
      break;
    }
  }

  if ((status == ConicStatus::NumericalError || status == ConicStatus::IterationLimit) && saved.merit < kInf) {
    x = saved.x, y = saved.y, z = saved.z, s = saved.s, tau = saved.tau, kappa = saved.kappa;
    pres = saved.pres, dres = saved.dres, gap = saved.gap;
    log.log(2, "conic   %s: the best iterate (reduced accuracy %.1e)\n", to_string(status), saved.merit);
    status = ConicStatus::Optimal;
    res.reduced_accuracy = true;
  }
  // ---- Unscale: x = D x~, y = E y~, z = F z~, s = s~ / F; the point divided by tau (a certificate: normalized).
  res.status = status;
  res.iterations = iter;
  const bool cert_p = status == ConicStatus::PrimalInfeasible, cert_d = status == ConicStatus::DualInfeasible;
  const double sc = cert_p ? -1.0 / (dot(b, y) + dot(h, z)) : cert_d ? -1.0 / dot(c, x) : 1.0 / tau;
  res.x.resize(static_cast<size_t>(n));
  res.y.resize(static_cast<size_t>(p));
  res.z.resize(static_cast<size_t>(m));
  res.s.resize(static_cast<size_t>(m));
  for (Index j = 0; j < n; ++j) res.x[j] = D[j] * x[j] * sc;
  for (Index i = 0; i < p; ++i) res.y[i] = E[i] * y[i] * sc;
  for (Index i = 0; i < m; ++i) res.z[i] = F[i] * z[i] * sc, res.s[i] = s[i] / F[i] * sc;
  double xpx_u = 0.0;
  if (hasP && !cert_p && !cert_d) {
    std::vector<double> pxu(static_cast<size_t>(n), 0.0);
    sym_mul_add(prob.P, res.x, pxu);
    xpx_u = dot(res.x, pxu);
  }
  res.pobj = 0.5 * xpx_u + dot(prob.c, res.x) + prob.offset;
  res.dobj = -0.5 * xpx_u - (dot(prob.b, res.y) + dot(prob.h, res.z)) + prob.offset;
  res.pres = pres, res.dres = dres, res.gap = gap;
  res.time = timer.elapsed();
  log.log(1, "Conic barrier: %s%s after %d iterations, objective %.10g (dual %.10g), %.2fs\n", to_string(status),
          res.reduced_accuracy ? " (reduced accuracy)" : "", iter, res.pobj, res.dobj, res.time);
  return res;
}

}  // namespace milp
