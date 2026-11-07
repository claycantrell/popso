#include "milp/solver/square_terms.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace milp {

namespace {

// Offsets beyond this are not used for the coarse grid (the rows' coefficients and right-hand sides grow
// with k^2); the tangents there are still valid, only weaker, and refine() adds windows where needed.
constexpr long long kGridMax = 1LL << 16;
constexpr long long kBoundCap = 1LL << 40;  // |x| beyond this counts as unbounded

long long clamp_ll(double v) {
  if (!(v > -static_cast<double>(kBoundCap))) return -kBoundCap;
  if (!(v < static_cast<double>(kBoundCap))) return kBoundCap;
  return std::llround(v);
}

}  // namespace

namespace {

bool integral(double v) { return std::isfinite(v) && v == std::round(v); }

// A continuous column whose value is an integer at every feasible point: it has coefficient +-1 in an equality
// row with an integer right-hand side whose other entries are integer columns with integer coefficients (an
// error e = sum x - T over integer production x). R: the matrix by rows.
bool implied_integer(const Model& m, const SparseMatrix& R, Index j) {
  for (Index k = m.A.begin(j); k < m.A.end(j); ++k) {
    const Index i = m.A.index[k];
    if (std::abs(m.A.value[k]) != 1.0 || m.row_lb[i] != m.row_ub[i] || !integral(m.row_lb[i])) continue;
    bool ok = true;
    for (Index t = R.begin(i); t < R.end(i) && ok; ++t) {
      const Index c = R.index[t];
      ok = c == j || (m.col_type[c] != VarType::Continuous && integral(R.value[t]));
    }
    if (ok) return true;
  }
  return false;
}

}  // namespace

bool SquareTerms::init(const Model& m, int window) {
  m_ = &m;
  window_ = std::max(1, window);
  terms_.clear();
  implied_.clear();
  if (m.has_logical() || !m.has_quadratic()) return false;
  const SparseMatrix Q = m.q_lower();
  const double sense = m.sense == ObjSense::Maximize ? -1.0 : 1.0;
  SparseMatrix R;  // by rows, built for the first continuous column
  for (Index j = 0; j < Q.num_major; ++j) {
    for (Index k = Q.begin(j); k < Q.end(j); ++k) {
      if (Q.index[k] != j) return false;  // a cross term: not separable
      const double q = Q.value[k];
      if (q == 0.0) continue;
      if (!(sense * q > 0.0)) return false;  // not convex in the model's sense
      const VarType t = m.col_type[j];
      if (t == VarType::Continuous) {
        if (R.num_major == 0) R = m.A.transpose();
        if (!implied_integer(m, R, j)) return false;
        implied_.push_back(j);
      } else if (t != VarType::Integer && t != VarType::Binary) {
        return false;
      }
      Term term;
      term.col = j;
      term.half_q = 0.5 * q;
      term.binary = t == VarType::Binary || (m.col_lb[j] >= 0.0 && m.col_ub[j] <= 1.0);
      terms_.push_back(term);
    }
  }
  if (terms_.empty()) return false;
  center_at(std::vector<double>(static_cast<size_t>(m.num_cols()), 0.0));
  return true;
}

void SquareTerms::add_window(Term& t, long long d) {
  // Tangent k is exact at d = k and d = k + 1: the offsets lo - 1 .. hi cover every d in [lo, hi].
  const long long a = std::max(t.lo - 1, d - window_), b = std::min(t.hi, d + window_);
  for (long long k = a; k <= b; ++k) t.ks.insert(k);
}

void SquareTerms::center_at(const std::vector<double>& x) {
  for (Term& t : terms_) {
    if (t.binary) continue;
    const long long xl = clamp_ll(std::ceil(m_->col_lb[t.col] - 1e-9));
    const long long xu = clamp_ll(std::floor(m_->col_ub[t.col] + 1e-9));
    const double v = x.empty() ? 0.0 : x[t.col];
    t.center = std::clamp(clamp_ll(std::isfinite(v) ? v : 0.0), xl, std::max(xl, xu));
    t.lo = xl == -kBoundCap ? -kBoundCap : xl - t.center;
    t.hi = xu == kBoundCap ? kBoundCap : xu - t.center;
    t.ks.clear();
    if (t.hi - t.lo <= 2LL * window_ + 2) {
      for (long long k = t.lo - 1; k <= t.hi; ++k) t.ks.insert(k);  // a narrow domain: every tangent
      continue;
    }
    add_window(t, 0);
    // A coarse grid beyond the window keeps the relaxation reasonable far from the center.
    for (long long g = 1; g <= kGridMax; g *= 2) {
      if (g <= window_) continue;
      for (long long k : {g, -g}) {
        if (k >= t.lo - 1 && k <= t.hi) t.ks.insert(k);
      }
    }
  }
}

namespace {

bool term_exact(long long d, const std::set<long long>& ks) {
  return d == 0 || ks.count(d) > 0 || ks.count(d - 1) > 0;  // t >= 0 makes d = 0 exact on its own
}

}  // namespace

bool SquareTerms::exact_at(const std::vector<double>& x) const {
  for (const Term& t : terms_) {
    if (!t.binary && !term_exact(std::llround(x[t.col]) - t.center, t.ks)) return false;
  }
  return true;
}

int SquareTerms::refine(const std::vector<double>& x) {
  int missed = 0;
  for (Term& t : terms_) {
    if (t.binary) continue;
    const long long d = std::llround(x[t.col]) - t.center;
    if (term_exact(d, t.ks)) continue;
    add_window(t, d);
    ++missed;
  }
  return missed;
}

size_t SquareTerms::num_general() const {
  size_t n = 0;
  for (const Term& t : terms_) n += !t.binary;
  return n;
}

size_t SquareTerms::num_tangents() const {
  size_t n = 0;
  for (const Term& t : terms_) n += t.ks.size();
  return n;
}

Model SquareTerms::build() const {
  const Model& m = *m_;
  Model out = m;
  for (Index j : implied_) out.col_type[j] = VarType::Integer;  // implied integer (init)
  out.q_row.clear();
  out.q_col.clear();
  out.q_val.clear();
  const Index n = m.num_cols(), m0 = m.num_rows();
  // New entries of the model's columns (tangent rows), by column; then the t columns.
  std::vector<std::vector<std::pair<Index, double>>> extra(static_cast<size_t>(n));
  std::vector<std::vector<Index>> t_rows;
  Index row = m0;
  for (const Term& t : terms_) {
    if (t.binary) {
      out.obj[t.col] += t.half_q;  // x^2 = x
      continue;
    }
    const double c = static_cast<double>(t.center);
    out.obj[t.col] += 2.0 * t.half_q * c;
    out.obj_offset -= t.half_q * c * c;
    out.col_names.push_back("popso_sq_t[" + m.col_names[t.col] + "]");
    out.obj.push_back(t.half_q);
    out.col_lb.push_back(0.0);
    out.col_ub.push_back(kInf);
    out.col_type.push_back(VarType::Continuous);
    t_rows.emplace_back();
    for (long long k : t.ks) {
      // t - (2k + 1) x >= -(2k + 1) c - k (k + 1), i.e. t >= (2k + 1) d - k (k + 1) with d = x - c.
      const double kk = static_cast<double>(k);
      out.row_names.push_back("popso_sq_tan[" + m.col_names[t.col] + "," + std::to_string(k) + "]");
      out.row_lb.push_back(-(2.0 * kk + 1.0) * c - kk * (kk + 1.0));
      out.row_ub.push_back(kInf);
      extra[t.col].emplace_back(row, -(2.0 * kk + 1.0));
      t_rows.back().push_back(row);
      ++row;
    }
  }
  SparseMatrix A;
  A.num_major = static_cast<Index>(out.obj.size());
  A.num_minor = row;
  A.start.assign(1, 0);
  for (Index j = 0; j < n; ++j) {
    for (Index k = m.A.begin(j); k < m.A.end(j); ++k) A.index.push_back(m.A.index[k]), A.value.push_back(m.A.value[k]);
    for (const auto& [r, v] : extra[j]) A.index.push_back(r), A.value.push_back(v);
    A.start.push_back(A.nnz());
  }
  for (const std::vector<Index>& rows : t_rows) {
    for (Index r : rows) A.index.push_back(r), A.value.push_back(1.0);
    A.start.push_back(A.nnz());
  }
  out.A = std::move(A);
  return out;
}

std::vector<double> SquareTerms::extend(const std::vector<double>& x) const {
  std::vector<double> y = x;
  for (const Term& t : terms_) {
    if (t.binary) continue;
    const double d = static_cast<double>(std::llround(x[t.col]) - t.center);
    y.push_back(d * d);
  }
  return y;
}

}  // namespace milp
