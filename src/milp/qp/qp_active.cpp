#include "milp/qp/qp_active.hpp"

#include <algorithm>
#include <cmath>

namespace milp {

SafeBound qp_safe_bound(const MipProblem& p, const std::vector<double>& qx, const std::vector<double>& x,
                        const std::vector<double>& y_in, const std::vector<double>& lb, const std::vector<double>& ub) {
  SafeBound r;
  const Index n = p.num_cols(), m = p.num_rows();
  std::vector<double> y = y_in;
  y.resize(static_cast<size_t>(m), 0.0);
  for (Index i = 0; i < m; ++i) {
    if (y[i] > 0.0 && !is_finite(p.row_lb[i])) y[i] = 0.0;
    if (y[i] < 0.0 && !is_finite(p.row_ub[i])) y[i] = 0.0;
  }
  r.d.assign(static_cast<size_t>(n), 0.0);
  double gmax = 1.0;
  for (Index j = 0; j < n; ++j) {
    double dj = p.cost[j] + qx[j];
    gmax = std::max(gmax, std::abs(dj));
    for (Index k = p.A.begin(j); k < p.A.end(j); ++k) dj -= p.A.value[k] * y[p.A.index[k]];
    r.d[j] = dj;
  }
  double xqx = 0.0;
  for (Index j = 0; j < n; ++j) xqx += x[j] * qx[j];
  double b = p.offset - 0.5 * xqx;
  double mag = std::abs(p.offset) + 0.5 * std::abs(xqx);
  for (Index i = 0; i < m; ++i) {
    const double t = y[i] > 0.0 ? y[i] * p.row_lb[i] : y[i] < 0.0 ? y[i] * p.row_ub[i] : 0.0;
    b += t;
    mag += std::abs(t);
  }
  // A reduced cost of round-off size on a column without the bound it points to is taken as 0
  // (the column is then valued at x_j): free columns have d_j = 0 only up to rounding.
  const double dtol = 1e-9 * gmax;
  for (Index j = 0; j < n; ++j) {
    const double dj = r.d[j];
    double t;
    if (dj > 0.0) t = is_finite(lb[j]) ? dj * lb[j] : dj <= dtol ? dj * x[j] : -kInf;
    else if (dj < 0.0) t = is_finite(ub[j]) ? dj * ub[j] : dj >= -dtol ? dj * x[j] : -kInf;
    else t = 0.0;
    if (!is_finite(t)) return r;  // bound -inf
    b += t;
    mag += std::abs(t);
  }
  r.bound = b - 1e-11 * mag;
  return r;
}

namespace {

// In-place Cholesky of the leading n x n block of a row-major matrix with leading dimension ld
// (lower triangle). Returns the first row whose pivot is <= rel_tol * its original diagonal (or
// not positive), or -1.
Index cholesky(double* a, Index n, Index ld, double rel_tol) {
  for (Index i = 0; i < n; ++i) {
    double* ai = a + static_cast<size_t>(i) * ld;
    for (Index j = 0; j < i; ++j) {
      const double* aj = a + static_cast<size_t>(j) * ld;
      double s = ai[j];
      for (Index k = 0; k < j; ++k) s -= ai[k] * aj[k];
      ai[j] = s / aj[j];
    }
    const double orig = ai[i];
    double s = orig;
    for (Index k = 0; k < i; ++k) s -= ai[k] * ai[k];
    if (!(s > rel_tol * std::abs(orig)) || !(s > 0.0)) return i;
    ai[i] = std::sqrt(s);
  }
  return -1;
}

void forward(const double* L, Index n, Index ld, double* b) {  // L b' = b
  for (Index i = 0; i < n; ++i) {
    const double* li = L + static_cast<size_t>(i) * ld;
    double s = b[i];
    for (Index k = 0; k < i; ++k) s -= li[k] * b[k];
    b[i] = s / li[i];
  }
}

void backward(const double* L, Index n, Index ld, double* b) {  // L^T b' = b
  for (Index i = n - 1; i >= 0; --i) {
    const double v = b[i] / L[static_cast<size_t>(i) * ld + i];
    b[i] = v;
    const double* li = L + static_cast<size_t>(i) * ld;
    for (Index k = 0; k < i; ++k) b[k] -= li[k] * v;
  }
}

}  // namespace

namespace qp_active_detail {
struct Border {
  bool add;  // added constraint (else: a removed base constraint)
  bool col;
  Index idx;
  int8_t st;              // added: its state (side)
  std::vector<double> w;  // K0^-1 b
};
}  // namespace qp_active_detail

// A factorization with its borders and the final working set of a solve (for warm starts).
struct ActiveSetQp::State {
  uint64_t id = 0;
  std::vector<int8_t> cst, rst, cst0, rst0;
  std::vector<Index> F, N0, R, posF, posN, posR;
  std::vector<double> L, V, S;
  std::vector<qp_active_detail::Border> borders;
  std::vector<double> C, Cf;
  std::vector<Index> Cpiv;
  bool c_ok = true;
  Index nf = 0, nn = 0, w = 0, len = 0;
};

ActiveSetQp::~ActiveSetQp() = default;

ActiveSetQp::ActiveSetQp(const MipProblem& p, const std::vector<double>& diag_shift) : p_(p) {
  n_ = p.num_cols();
  m_ = p.num_rows();
  Qd_.assign(static_cast<size_t>(n_) * n_, 0.0);
  for (Index j = 0; j < p.Q.num_major; ++j) {
    for (Index k = p.Q.begin(j); k < p.Q.end(j); ++k) {
      const Index i = p.Q.index[k];
      Qd_[static_cast<size_t>(i) * n_ + j] += p.Q.value[k];
      if (i != j) Qd_[static_cast<size_t>(j) * n_ + i] += p.Q.value[k];
    }
  }
  if (!diag_shift.empty())
    for (Index j = 0; j < n_; ++j) Qd_[static_cast<size_t>(j) * n_ + j] += diag_shift[j];
  double sum = 0.0;
  Index cnt = 0;
  for (Index j = 0; j < n_; ++j) {
    const double q = Qd_[static_cast<size_t>(j) * n_ + j];
    if (q > 0.0) sum += q, ++cnt;
  }
  const double qavg = cnt > 0 ? sum / cnt : 1.0;
  eps_.resize(static_cast<size_t>(n_));
  for (Index j = 0; j < n_; ++j) eps_[j] = 1e-8 * (std::max(0.0, Qd_[static_cast<size_t>(j) * n_ + j]) + qavg);
  rows_ = p.A.transpose();
  row_norm_.assign(static_cast<size_t>(m_), 0.0);
  for (Index i = 0; i < m_; ++i) {
    double s = 0.0;
    for (Index k = rows_.begin(i); k < rows_.end(i); ++k) s += rows_.value[k] * rows_.value[k];
    row_norm_[i] = std::max(1e-300, std::sqrt(s));
  }
}

void ActiveSetQp::q_times(const std::vector<double>& x, std::vector<double>& out) const {
  out.assign(static_cast<size_t>(n_), 0.0);
  for (Index i = 0; i < n_; ++i) {
    const double* qi = Qd_.data() + static_cast<size_t>(i) * n_;
    double s = 0.0;
    for (Index j = 0; j < n_; ++j) s += qi[j] * x[j];
    out[i] = s;
  }
}

WorkingSet ActiveSetQp::guess(const std::vector<double>& x, const std::vector<double>& lb,
                              const std::vector<double>& ub, double tol) const {
  WorkingSet ws;
  for (Index j = 0; j < n_; ++j) {
    if (is_finite(lb[j]) && x[j] - lb[j] <= tol * (1.0 + std::abs(lb[j]))) ws.push_back(4 * j);
    else if (is_finite(ub[j]) && ub[j] - x[j] <= tol * (1.0 + std::abs(ub[j]))) ws.push_back(4 * j + 1);
  }
  for (Index i = 0; i < m_; ++i) {
    double act = 0.0;
    for (Index k = rows_.begin(i); k < rows_.end(i); ++k) act += rows_.value[k] * x[rows_.index[k]];
    const double lo = p_.row_lb[i], hi = p_.row_ub[i];
    if (is_finite(lo) && act - lo <= tol * (1.0 + std::abs(lo)) * row_norm_[i]) ws.push_back(4 * i + 2);
    else if (is_finite(hi) && hi - act <= tol * (1.0 + std::abs(hi)) * row_norm_[i]) ws.push_back(4 * i + 3);
  }
  return ws;
}

namespace {

// Dense LU with partial pivoting of a k x k row-major matrix (in place); false if (nearly) singular.
bool lu_factor(std::vector<double>& a, std::vector<Index>& piv, Index k) {
  piv.resize(static_cast<size_t>(k));
  double amax = 0.0;
  for (const double v : a) amax = std::max(amax, std::abs(v));
  for (Index c = 0; c < k; ++c) {
    Index p = c;
    for (Index r = c + 1; r < k; ++r)
      if (std::abs(a[static_cast<size_t>(r) * k + c]) > std::abs(a[static_cast<size_t>(p) * k + c])) p = r;
    piv[c] = p;
    if (!(std::abs(a[static_cast<size_t>(p) * k + c]) > 1e-13 * amax)) return false;
    if (p != c)
      for (Index t = 0; t < k; ++t) std::swap(a[static_cast<size_t>(c) * k + t], a[static_cast<size_t>(p) * k + t]);
    const double d = a[static_cast<size_t>(c) * k + c];
    for (Index r = c + 1; r < k; ++r) {
      const double f = a[static_cast<size_t>(r) * k + c] / d;
      a[static_cast<size_t>(r) * k + c] = f;
      if (f != 0.0)
        for (Index t = c + 1; t < k; ++t) a[static_cast<size_t>(r) * k + t] -= f * a[static_cast<size_t>(c) * k + t];
    }
  }
  return true;
}

void lu_solve(const std::vector<double>& a, const std::vector<Index>& piv, Index k, std::vector<double>& b) {
  // The factorization swapped whole rows (multipliers included): permute first, then eliminate.
  for (Index c = 0; c < k; ++c) std::swap(b[c], b[piv[c]]);
  for (Index c = 0; c < k; ++c)
    for (Index r = c + 1; r < k; ++r) b[r] -= a[static_cast<size_t>(r) * k + c] * b[c];
  for (Index c = k - 1; c >= 0; --c) {
    double s = b[c];
    for (Index t = c + 1; t < k; ++t) s -= a[static_cast<size_t>(c) * k + t] * b[t];
    b[c] = s / a[static_cast<size_t>(c) * k + c];
  }
}

}  // namespace

// The KKT system of a working set W (columns held at a bound N, active rows R), with H = Q +
// diag(eps) and multipliers nu (signed: d_j for bounds, y_i for rows),
//
//   [ H     -N ] [x ]   [r]
//   [ -N^T   0 ] [nu] = [s],
//
// is factored once for a base working set W0 (bounds eliminated: Cholesky of H_FF, Schur
// complement S = A_RF H_FF^-1 A_RF^T); later changes of W are borders of that system (an added
// constraint c: column [-n_c; 0]; a removed base constraint k: column [0; e_k], which forces
// nu_k = 0 and frees its equation), solved through their dense Schur complement C. After
// kMaxBorders changes the base is refactored.
ActiveSetResult ActiveSetQp::solve(const std::vector<double>& lb, const std::vector<double>& ub,
                                   const WorkingSet& start, const std::vector<double>& center_in,
                                   long long flop_budget, uint64_t start_state) {
  using St = ActiveSetResult::Status;
  constexpr Index kMaxBorders = 40;
  ActiveSetResult res;
  const Index n = n_, m = m_;
  const std::vector<double>& rlb = p_.row_lb;
  const std::vector<double>& rub = p_.row_ub;
  const SparseMatrix& A = p_.A;
  // Column states: 0 free, 1 at lower, 2 at upper, 3 fixed (lb == ub). Rows: 0 inactive,
  // 1 lower side, 2 upper side, 3 equality.
  std::vector<int8_t> cst(static_cast<size_t>(n), 0), rst(static_cast<size_t>(m), 0);
  for (Index j = 0; j < n; ++j) {
    if (lb[j] > ub[j]) {
      res.status = St::Infeasible;
      return res;
    }
    if (lb[j] == ub[j]) cst[j] = 3;
  }
  for (const int32_t e : start) {
    const Index k = e >> 2;
    const int t = e & 3;
    if (t < 2) {
      if (k >= n || cst[k] == 3) continue;
      if (t == 0 && is_finite(lb[k])) cst[k] = 1;
      else if (t == 1 && is_finite(ub[k])) cst[k] = 2;
    } else {
      if (k >= m) continue;
      if (is_finite(rlb[k]) && rlb[k] == rub[k]) rst[k] = 3;
      else if (t == 2 && is_finite(rlb[k])) rst[k] = 1;
      else if (t == 3 && is_finite(rub[k])) rst[k] = 2;
    }
  }
  std::vector<double> center = center_in;
  center.resize(static_cast<size_t>(n), 0.0);
  for (Index j = 0; j < n; ++j) center[j] = std::min(std::max(center[j], lb[j]), ub[j]);
  double hmax = 0.0, gscale = 1.0;
  for (Index j = 0; j < n; ++j) {
    hmax = std::max(hmax, Qd_[static_cast<size_t>(j) * n + j] + eps_[j]);
    gscale = std::max(gscale, std::abs(p_.cost[j]));
  }
  auto col_value = [&](Index j, int8_t st) { return st == 2 ? ub[j] : lb[j]; };
  auto row_value = [&](Index i, int8_t st) { return st == 2 ? rub[i] : rlb[i]; };
  // The same constraint (held at the same side; a fixed column or an equality row matches either).
  auto same = [](int8_t a, int8_t b) { return a != 0 && b != 0 && (a == b || a == 3 || b == 3); };

  // ---- Base factorization.
  std::vector<int8_t> cst0, rst0;
  std::vector<Index> F, N0, R, posF(static_cast<size_t>(n), -1), posN(static_cast<size_t>(n), -1),
      posR(static_cast<size_t>(m), -1);
  std::vector<double> L, V, S, Lc;
  std::vector<Index> keep;
  Index nf = 0, nn = 0, w = 0, len = 0;
  long long flops = 0;
  // Borders.
  using Border = qp_active_detail::Border;
  std::vector<Border> borders;
  std::vector<double> C, Cf;
  std::vector<Index> Cpiv;
  bool c_ok = true;

  auto factor = [&](bool allow_drop) -> bool {
    borders.clear();
    F.clear();
    N0.clear();
    R.clear();
    for (Index j = 0; j < n; ++j) {
      posF[j] = posN[j] = -1;
      if (cst[j] == 0) posF[j] = static_cast<Index>(F.size()), F.push_back(j);
      else posN[j] = static_cast<Index>(N0.size()), N0.push_back(j);
    }
    for (Index i = 0; i < m; ++i) posR[i] = -1;
    for (Index i = 0; i < m; ++i)
      if (rst[i] != 0) R.push_back(i);
    nf = static_cast<Index>(F.size());
    nn = static_cast<Index>(N0.size());
    w = static_cast<Index>(R.size());
    L.assign(static_cast<size_t>(nf) * nf, 0.0);
    for (Index a = 0; a < nf; ++a) {
      const double* qa = Qd_.data() + static_cast<size_t>(F[a]) * n;
      double* la = L.data() + static_cast<size_t>(a) * nf;
      for (Index b = 0; b <= a; ++b) la[b] = qa[F[b]];
      la[a] += eps_[F[a]];
    }
    flops += static_cast<long long>(nf) * nf * nf / 3 + nf;
    if (cholesky(L.data(), nf, nf, 1e-14) >= 0) return false;
    V.assign(static_cast<size_t>(w) * nf, 0.0);
    for (Index r = 0; r < w; ++r) {
      double* vr = V.data() + static_cast<size_t>(r) * nf;
      const Index i = R[r];
      for (Index k = rows_.begin(i); k < rows_.end(i); ++k) {
        const Index j = rows_.index[k];
        if (posF[j] >= 0) vr[posF[j]] = rows_.value[k];
      }
      forward(L.data(), nf, nf, vr);
    }
    S.assign(static_cast<size_t>(w) * w, 0.0);
    for (Index r = 0; r < w; ++r) {
      const double* vr = V.data() + static_cast<size_t>(r) * nf;
      for (Index s = 0; s <= r; ++s) {
        const double* vs = V.data() + static_cast<size_t>(s) * nf;
        double t = 0.0;
        for (Index k = 0; k < nf; ++k) t += vr[k] * vs[k];
        S[static_cast<size_t>(r) * w + s] = t;
      }
    }
    flops += static_cast<long long>(w) * nf * nf + static_cast<long long>(w) * w * nf / 2 +
             static_cast<long long>(w) * w * w / 3;
    // Cholesky of S row by row; a row dependent on the earlier ones (pivot ~ 0) leaves the
    // working set when allowed.
    Lc.assign(static_cast<size_t>(w) * w, 0.0);
    keep.clear();
    for (Index r = 0; r < w; ++r) {
      const double* sr = S.data() + static_cast<size_t>(r) * w;
      double* lr = Lc.data() + keep.size() * static_cast<size_t>(w);
      const Index t_end = static_cast<Index>(keep.size());
      double piv = sr[r];
      for (Index t = 0; t < t_end; ++t) {
        const double* lt = Lc.data() + static_cast<size_t>(t) * w;
        double s = sr[keep[t]];
        for (Index k = 0; k < t; ++k) s -= lr[k] * lt[k];
        lr[t] = s / lt[t];
        piv -= lr[t] * lr[t];
      }
      if (!(piv > 1e-11 * std::abs(sr[r])) || !(piv > 0.0)) {
        if (!allow_drop) return false;
        rst[R[r]] = 0;
        std::fill(lr, lr + w, 0.0);
        continue;
      }
      lr[t_end] = std::sqrt(piv);
      keep.push_back(r);
    }
    const Index w2 = static_cast<Index>(keep.size());
    for (Index t = 0; t < w2; ++t) {
      R[t] = R[keep[t]];
      if (t != keep[t])
        std::copy_n(V.data() + static_cast<size_t>(keep[t]) * nf, nf, V.data() + static_cast<size_t>(t) * nf);
    }
    R.resize(static_cast<size_t>(w2));
    V.resize(static_cast<size_t>(w2) * nf);
    S.assign(static_cast<size_t>(w2) * w2, 0.0);
    for (Index t = 0; t < w2; ++t)
      std::copy_n(Lc.data() + static_cast<size_t>(t) * w, t + 1, S.data() + static_cast<size_t>(t) * w2);
    w = w2;
    for (Index r = 0; r < w; ++r) posR[R[r]] = r;
    cst0 = cst;
    rst0 = rst;
    len = n + nn + w;
    C.clear();
    c_ok = true;
    return true;
  };

  // K0 solve. Layout of r / out: [x (n) | bound multipliers (nn) | row multipliers (w)].
  std::vector<double> s_rf, s_y, s_yfull(static_cast<size_t>(m), 0.0);
  std::vector<Index> s_nz;
  auto solve0 = [&](const std::vector<double>& r, std::vector<double>& out) {
    out.assign(static_cast<size_t>(len), 0.0);
    s_nz.clear();
    for (Index t = 0; t < nn; ++t) {
      const double v = -r[n + t];
      out[N0[t]] = v;
      if (v != 0.0) s_nz.push_back(N0[t]);
    }
    s_rf.assign(static_cast<size_t>(nf), 0.0);
    for (Index a = 0; a < nf; ++a) s_rf[a] = r[F[a]];
    for (const Index j : s_nz) {
      const double xj = out[j];
      const double* qj = Qd_.data() + static_cast<size_t>(j) * n;
      for (Index a = 0; a < nf; ++a) s_rf[a] -= qj[F[a]] * xj;
    }
    forward(L.data(), nf, nf, s_rf.data());
    s_y.assign(static_cast<size_t>(w), 0.0);
    for (Index rr = 0; rr < w; ++rr) {
      const Index i = R[rr];
      double b = -r[n + nn + rr];
      if (!s_nz.empty())
        for (Index k = rows_.begin(i); k < rows_.end(i); ++k)
          if (posN[rows_.index[k]] >= 0) b -= rows_.value[k] * out[rows_.index[k]];
      const double* vr = V.data() + static_cast<size_t>(rr) * nf;
      for (Index k = 0; k < nf; ++k) b -= vr[k] * s_rf[k];
      s_y[rr] = b;
    }
    forward(S.data(), w, w, s_y.data());
    backward(S.data(), w, w, s_y.data());
    for (Index rr = 0; rr < w; ++rr) {
      const double* vr = V.data() + static_cast<size_t>(rr) * nf;
      const double yr = s_y[rr];
      if (yr != 0.0)
        for (Index k = 0; k < nf; ++k) s_rf[k] += yr * vr[k];
    }
    backward(L.data(), nf, nf, s_rf.data());
    for (Index a = 0; a < nf; ++a) out[F[a]] = s_rf[a];
    for (Index rr = 0; rr < w; ++rr) s_yfull[R[rr]] = s_y[rr], out[n + nn + rr] = s_y[rr];
    for (Index t = 0; t < nn; ++t) {
      const Index j = N0[t];
      const double* qj = Qd_.data() + static_cast<size_t>(j) * n;
      double hx = eps_[j] * out[j];
      for (Index a = 0; a < nf; ++a) hx += qj[F[a]] * s_rf[a];
      for (const Index k : s_nz) hx += qj[k] * out[k];
      for (Index k = A.begin(j); k < A.end(j); ++k) hx -= A.value[k] * s_yfull[A.index[k]];
      out[n + t] = hx - r[j];
    }
    for (Index rr = 0; rr < w; ++rr) s_yfull[R[rr]] = 0.0;
    flops += 2LL * nf * nf + 2LL * w * nf + static_cast<long long>(w) * w + static_cast<long long>(nn) * nf +
             static_cast<long long>(s_nz.size()) * (nf + nn) + 2LL * A.nnz();
  };
  // b^T u and the border column b.
  auto bdot = [&](const Border& b, const std::vector<double>& u) -> double {
    if (b.add) {
      if (b.col) return -u[b.idx];
      double s = 0.0;
      for (Index k = rows_.begin(b.idx); k < rows_.end(b.idx); ++k) s -= rows_.value[k] * u[rows_.index[k]];
      return s;
    }
    return b.col ? u[n + posN[b.idx]] : u[n + nn + posR[b.idx]];
  };
  std::vector<double> bvec;
  auto set_border = [&](Border& b) {
    bvec.assign(static_cast<size_t>(len), 0.0);
    if (b.add) {
      if (b.col) bvec[b.idx] = -1.0;
      else
        for (Index k = rows_.begin(b.idx); k < rows_.end(b.idx); ++k) bvec[rows_.index[k]] = -rows_.value[k];
    } else {
      bvec[b.col ? n + posN[b.idx] : n + nn + posR[b.idx]] = 1.0;
    }
    solve0(bvec, b.w);
  };
  auto rebuild_c = [&]() {
    const Index k = static_cast<Index>(borders.size());
    C.assign(static_cast<size_t>(k) * k, 0.0);
    for (Index a = 0; a < k; ++a)
      for (Index b = 0; b < k; ++b) C[static_cast<size_t>(a) * k + b] = bdot(borders[a], borders[b].w);
    Cf = C;
    c_ok = k == 0 || lu_factor(Cf, Cpiv, k);
    flops += static_cast<long long>(k) * k * (k + 10);
  };
  // Bordered solve: rhs r (K0 layout) and s (one entry per border). out: K0 layout, v: borders.
  std::vector<double> u0, tv;
  auto bsolve = [&](const std::vector<double>& r, const std::vector<double>& s, std::vector<double>& out,
                    std::vector<double>& v) {
    solve0(r, out);
    const Index k = static_cast<Index>(borders.size());
    v.assign(static_cast<size_t>(k), 0.0);
    if (k == 0) return;
    for (Index a = 0; a < k; ++a) v[a] = bdot(borders[a], out) - s[a];
    lu_solve(Cf, Cpiv, k, v);
    for (Index a = 0; a < k; ++a) {
      const double va = v[a];
      const std::vector<double>& wa = borders[a].w;
      for (Index t = 0; t < len; ++t) out[t] -= va * wa[t];
    }
    flops += static_cast<long long>(k) * (2 * len + k);
  };
  // Multipliers of the current working set from a bordered solution.
  auto unpack = [&](const std::vector<double>& out, const std::vector<double>& v, std::vector<double>& xo,
                    std::vector<double>& mcol, std::vector<double>& mrow) {
    xo.assign(out.begin(), out.begin() + n);
    mcol.assign(static_cast<size_t>(n), 0.0);
    mrow.assign(static_cast<size_t>(m), 0.0);
    for (Index t = 0; t < nn; ++t) {
      const Index j = N0[t];
      if (same(cst[j], cst0[j])) mcol[j] = out[n + t];
    }
    for (Index rr = 0; rr < w; ++rr) {
      const Index i = R[rr];
      if (same(rst[i], rst0[i])) mrow[i] = out[n + nn + rr];
    }
    for (size_t a = 0; a < borders.size(); ++a) {
      const Border& b = borders[a];
      if (!b.add) continue;
      if (b.col) mcol[b.idx] = v[a];
      else mrow[b.idx] = v[a];
    }
  };

  // Primal / dual state.
  std::vector<double> x(static_cast<size_t>(n), 0.0), y(static_cast<size_t>(m), 0.0), d(static_cast<size_t>(n), 0.0);
  std::vector<double> act(static_cast<size_t>(m), 0.0);
  // Pending constraint (being added): column or row index, orientation s (normal s e_j or s a_i),
  // and its current multiplier.
  bool pend = false, pend_col = false;
  Index pend_idx = -1;
  double pend_s = 1.0, pend_u = 0.0;
  std::vector<double> rhs, sv, out, vv;
  // The minimizer on the working set (with the pending constraint's multiplier).
  auto eqp = [&]() {
    rhs.assign(static_cast<size_t>(len), 0.0);
    for (Index j = 0; j < n; ++j) rhs[j] = -(p_.cost[j] - eps_[j] * center[j]);
    if (pend) {
      if (pend_col) rhs[pend_idx] += pend_u * pend_s;
      else
        for (Index k = rows_.begin(pend_idx); k < rows_.end(pend_idx); ++k)
          rhs[rows_.index[k]] += pend_u * pend_s * rows_.value[k];
    }
    for (Index t = 0; t < nn; ++t) rhs[n + t] = -col_value(N0[t], cst0[N0[t]]);
    for (Index rr = 0; rr < w; ++rr) rhs[n + nn + rr] = -row_value(R[rr], rst0[R[rr]]);
    sv.assign(borders.size(), 0.0);
    for (size_t a = 0; a < borders.size(); ++a) {
      const Border& b = borders[a];
      if (b.add) sv[a] = b.col ? -col_value(b.idx, b.st) : -row_value(b.idx, b.st);
    }
    bsolve(rhs, sv, out, vv);
    unpack(out, vv, x, d, y);
    for (Index j = 0; j < n; ++j)
      if (cst[j] != 0) x[j] = col_value(j, cst[j]);
    std::fill(act.begin(), act.end(), 0.0);
    for (Index j = 0; j < n; ++j)
      if (x[j] != 0.0)
        for (Index k = A.begin(j); k < A.end(j); ++k) act[A.index[k]] += A.value[k] * x[j];
    flops += 2LL * A.nnz() + n;
  };
  // Working set changes.
  auto find_border = [&](bool add, bool col, Index idx) -> Index {
    for (size_t a = 0; a < borders.size(); ++a)
      if (borders[a].add == add && borders[a].col == col && borders[a].idx == idx) return static_cast<Index>(a);
    return -1;
  };
  auto push_border = [&](Border b) {
    set_border(b);
    borders.push_back(std::move(b));
    rebuild_c();
  };
  auto erase_border = [&](Index a) {
    if (a < 0) {
      c_ok = false;  // inconsistent borders: refactor
      return;
    }
    borders.erase(borders.begin() + a);
    rebuild_c();
  };
  auto add_constraint = [&](bool col, Index idx, int8_t st) {
    // Re-adding a removed base constraint deletes its border.
    const Index a = same(col ? cst0[idx] : rst0[idx], st) ? find_border(false, col, idx) : -1;
    if (col) cst[idx] = st;
    else rst[idx] = st;
    if (a >= 0) erase_border(a);
    else push_border(Border{true, col, idx, st, {}});
  };
  auto remove_constraint = [&](bool col, Index idx) {
    if (col) {
      const bool base = same(cst0[idx], cst[idx]);
      cst[idx] = 0;
      if (base) push_border(Border{false, true, idx, 0, {}});
      else erase_border(find_border(true, true, idx));
    } else {
      const bool base = same(rst0[idx], rst[idx]);
      rst[idx] = 0;
      if (base) push_border(Border{false, false, idx, 0, {}});
      else erase_border(find_border(true, false, idx));
    }
  };
  const double dtol = 1e-9 * gscale;
  const double ptol = 1e-9;

  // Start. A cached factorization of the start working set is reused: the differences that
  // the new bounds make (columns fixed now, ...) become borders.
  bool loaded = false;
  for (auto it = cache_.begin(); it != cache_.end() && start_state != 0; ++it) {
    State& c = **it;
    if (c.id != start_state) continue;
    const std::vector<int8_t> cst_d = cst, rst_d = rst;
    cst = c.cst, rst = c.rst, cst0 = c.cst0, rst0 = c.rst0;
    F = c.F, N0 = c.N0, R = c.R, posF = c.posF, posN = c.posN, posR = c.posR;
    L = c.L, V = c.V, S = c.S, borders = c.borders, C = c.C, Cf = c.Cf, Cpiv = c.Cpiv;
    c_ok = c.c_ok, nf = c.nf, nn = c.nn, w = c.w, len = c.len;
    loaded = true;
    for (Index j = 0; j < n && c_ok && static_cast<Index>(borders.size()) <= kMaxBorders; ++j) {
      if (cst[j] == cst_d[j]) continue;
      if (same(cst[j], cst_d[j])) {
        // Same side; the state code changes (fixed <-> at a bound) but not the matrix.
        const bool base = same(cst0[j], cst[j]);
        cst[j] = cst_d[j];
        if (base) cst0[j] = cst_d[j];
        else
          for (Border& b : borders)
            if (b.add && b.col && b.idx == j) b.st = cst_d[j];
        continue;
      }
      if (cst[j] != 0) remove_constraint(true, j);
      if (cst_d[j] != 0) add_constraint(true, j, cst_d[j]);
    }
    for (Index i = 0; i < m && c_ok && static_cast<Index>(borders.size()) <= kMaxBorders; ++i) {
      if (rst[i] == rst_d[i] || same(rst[i], rst_d[i])) continue;
      if (rst[i] != 0) remove_constraint(false, i);
      if (rst_d[i] != 0) add_constraint(false, i, rst_d[i]);
    }
    if (!c_ok || static_cast<Index>(borders.size()) > kMaxBorders) {
      cst = cst_d, rst = rst_d;
      loaded = false;
    }
    // Most recently used last (a node's state stays cached while its children are tried).
    std::rotate(it, it + 1, cache_.end());
    break;
  }
  if (!loaded && !factor(true)) {
    res.reason = "start factorization";
    return res;
  }
  eqp();
  // A dual feasible working set: constraints whose multiplier has the wrong sign are dropped
  // (a few as borders, else all at once with a new factorization) until there are none.
  for (int guard = 0;; ++guard) {
    if (guard > n + m || flops > flop_budget) {
      res.reason = "start repair";
      res.flops = flops;
      return res;
    }
    std::vector<std::pair<bool, Index>> wrong;
    for (Index j = 0; j < n; ++j) {
      const double v = cst[j] == 1 ? -d[j] : cst[j] == 2 ? d[j] : 0.0;
      if (v > dtol) wrong.push_back({true, j});
    }
    for (Index i = 0; i < m; ++i) {
      const double v = rst[i] == 1 ? -y[i] : rst[i] == 2 ? y[i] : 0.0;
      if (v > dtol) wrong.push_back({false, i});
    }
    if (wrong.empty()) break;
    if (wrong.size() <= 4 && static_cast<Index>(borders.size() + wrong.size()) <= kMaxBorders) {
      for (const auto& [col, k] : wrong) remove_constraint(col, k);
      if (c_ok) {
        eqp();
        continue;
      }
    } else {
      for (const auto& [col, k] : wrong) (col ? cst[k] : rst[k]) = 0;
    }
    if (!factor(true)) {
      res.reason = "start factorization";
      return res;
    }
    eqp();
  }

  const int max_iter = 5 * (n + m) + 50;
  std::vector<double> z, dcol, drow, dout, dv;
  int prox_rounds = 0;
  while (true) {
    if (res.iterations > max_iter || flops > flop_budget) {
      res.reason = res.iterations > max_iter ? "iteration limit" : "work limit";
      res.flops = flops;
      return res;
    }
    if (!pend) {
      // Most violated constraint (scaled by the row norm).
      double best = 0.0;
      for (Index j = 0; j < n; ++j) {
        if (cst[j] != 0) continue;
        const double vl = lb[j] - x[j], vu = x[j] - ub[j];
        if (vl > ptol * (1.0 + std::abs(lb[j])) && vl > best) best = vl, pend_idx = j, pend_col = true, pend_s = 1.0;
        if (vu > ptol * (1.0 + std::abs(ub[j])) && vu > best) best = vu, pend_idx = j, pend_col = true, pend_s = -1.0;
      }
      for (Index i = 0; i < m; ++i) {
        if (rst[i] != 0) continue;
        const double vl = (rlb[i] - act[i]) / row_norm_[i], vu = (act[i] - rub[i]) / row_norm_[i];
        if (vl > ptol * (1.0 + std::abs(rlb[i])) && vl > best) best = vl, pend_idx = i, pend_col = false, pend_s = 1.0;
        if (vu > ptol * (1.0 + std::abs(rub[i])) && vu > best) best = vu, pend_idx = i, pend_col = false, pend_s = -1.0;
      }
      if (best <= 0.0) {
        // Optimal for the regularized QP. Safe bound of the true one; one more round around the
        // new center if the regularization costs a noticeable part of it.
        res.x = x;
        res.y = y;
        std::vector<double> qx;
        q_times(x, qx);
        res.safe = qp_safe_bound(p_, qx, x, y, lb, ub);
        double obj = p_.offset;
        for (Index j = 0; j < n; ++j) obj += x[j] * (p_.cost[j] + 0.5 * qx[j]);
        res.objective = obj;
        flops += 2LL * n * n;
        const double loss = obj - res.safe.bound;
        if (prox_rounds < 2 && !(loss <= 1e-9 * (1.0 + std::abs(obj)))) {
          ++prox_rounds;
          center = x;
          eqp();
          continue;
        }
        break;
      }
      pend = true;
      pend_u = 0.0;
    }
    ++res.iterations;
    // Direction: H z - N r = n_p, N^T z = 0.
    rhs.assign(static_cast<size_t>(len), 0.0);
    if (pend_col) rhs[pend_idx] = pend_s;
    else
      for (Index k = rows_.begin(pend_idx); k < rows_.end(pend_idx); ++k) rhs[rows_.index[k]] = pend_s * rows_.value[k];
    sv.assign(borders.size(), 0.0);
    bsolve(rhs, sv, dout, dv);
    unpack(dout, dv, z, dcol, drow);
    double rate = 0.0, npn = 0.0, rmax = 0.0;
    if (pend_col) {
      rate = pend_s * z[pend_idx];
      npn = 1.0;
    } else {
      for (Index k = rows_.begin(pend_idx); k < rows_.end(pend_idx); ++k) {
        rate += pend_s * rows_.value[k] * z[rows_.index[k]];
        npn += rows_.value[k] * rows_.value[k];
      }
    }
    for (Index j = 0; j < n; ++j) rmax = std::max(rmax, std::abs(dcol[j]));
    for (Index i = 0; i < m; ++i) rmax = std::max(rmax, std::abs(drow[i]));
    // Step lengths.
    const double viol = pend_col ? (pend_s > 0 ? lb[pend_idx] - x[pend_idx] : x[pend_idx] - ub[pend_idx])
                                 : (pend_s > 0 ? rlb[pend_idx] - act[pend_idx] : act[pend_idx] - rub[pend_idx]);
    double t2 = kInf;
    if (rate > 1e-13 * npn / std::max(hmax, 1e-300)) t2 = std::max(0.0, viol) / rate;
    double t1 = kInf;
    Index drop = -1;
    bool drop_is_col = false;
    const double rtol = 1e-11 * std::max(rmax, 1e-300);
    for (Index i = 0; i < m; ++i) {
      const double ri = drow[i];
      double t = kInf;
      if (rst[i] == 1 && ri < -rtol) t = std::max(0.0, y[i]) / -ri;
      else if (rst[i] == 2 && ri > rtol) t = std::max(0.0, -y[i]) / ri;
      if (t < t1) t1 = t, drop = i, drop_is_col = false;
    }
    for (Index j = 0; j < n; ++j) {
      double t = kInf;
      if (cst[j] == 1 && dcol[j] < -rtol) t = std::max(0.0, d[j]) / -dcol[j];
      else if (cst[j] == 2 && dcol[j] > rtol) t = std::max(0.0, -d[j]) / dcol[j];
      if (t < t1) t1 = t, drop = j, drop_is_col = true;
    }
    if (!is_finite(t1) && !is_finite(t2)) {
      // Infeasible: verify the Farkas certificate of the row multipliers (the active rows'
      // direction, the pending row with 1) against the bounds.
      std::vector<double> yr = drow;
      if (!pend_col) yr[pend_idx] += pend_s;
      double Y = 0.0, U = 0.0, mag = 0.0;
      bool ok = true;
      for (Index i = 0; i < m && ok; ++i) {
        if (yr[i] > 0.0) ok = is_finite(rlb[i]), Y += yr[i] * rlb[i], mag += std::abs(yr[i] * rlb[i]);
        else if (yr[i] < 0.0) ok = is_finite(rub[i]), Y += yr[i] * rub[i], mag += std::abs(yr[i] * rub[i]);
      }
      for (Index j = 0; j < n && ok; ++j) {
        double v = 0.0;
        for (Index k = A.begin(j); k < A.end(j); ++k) v += A.value[k] * yr[A.index[k]];
        if (v > 0.0) ok = is_finite(ub[j]), U += v * ub[j], mag += std::abs(v * ub[j]);
        else if (v < 0.0) ok = is_finite(lb[j]), U += v * lb[j], mag += std::abs(v * lb[j]);
      }
      res.flops = flops;
      if (ok && Y - U > 1e-9 * std::max(1.0, mag)) {
        res.status = St::Infeasible;
        return res;
      }
      res.reason = "unverified infeasibility";
      return res;
    }
    if (t2 <= t1) {
      const int8_t st = pend_col ? (pend_s > 0 ? 1 : 2) : (rlb[pend_idx] == rub[pend_idx] ? 3 : pend_s > 0 ? 1 : 2);
      pend = false;
      pend_u = 0.0;
      add_constraint(pend_col, pend_idx, st);
    } else {
      pend_u += t1;
      remove_constraint(drop_is_col, drop);
    }
    if (!c_ok || static_cast<Index>(borders.size()) > kMaxBorders) {
      if (!factor(false)) {
        res.reason = "dependent working set";
        res.flops = flops;
        return res;
      }
    }
    eqp();
  }
  res.status = St::Optimal;
  res.flops = flops;
  {
    auto c = std::make_unique<State>();
    c->id = next_state_++;
    c->cst = cst, c->rst = rst, c->cst0 = std::move(cst0), c->rst0 = std::move(rst0);
    c->F = std::move(F), c->N0 = std::move(N0), c->R = std::move(R);
    c->posF = std::move(posF), c->posN = std::move(posN), c->posR = std::move(posR);
    c->L = std::move(L), c->V = std::move(V), c->S = std::move(S), c->borders = std::move(borders);
    c->C = std::move(C), c->Cf = std::move(Cf), c->Cpiv = std::move(Cpiv);
    c->c_ok = c_ok, c->nf = nf, c->nn = nn, c->w = w, c->len = len;
    res.state = c->id;
    if (cache_.size() >= 4) cache_.erase(cache_.begin());
    cache_.push_back(std::move(c));
  }
  for (Index j = 0; j < n; ++j) {
    if (cst[j] == 1 || cst[j] == 3) res.ws.push_back(4 * j);
    else if (cst[j] == 2) res.ws.push_back(4 * j + 1);
  }
  for (Index i = 0; i < m; ++i) {
    if (rst[i] == 1 || rst[i] == 3) res.ws.push_back(4 * i + 2);
    else if (rst[i] == 2) res.ws.push_back(4 * i + 3);
  }
  return res;
}

}  // namespace milp
