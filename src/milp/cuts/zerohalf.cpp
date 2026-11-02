// Zero-half cuts ({0, 1/2}-Chvatal-Gomory cuts).
//
// Every problem row side whose integer coefficients (and right-hand side, after
// replacing the continuous columns by the bound that relaxes the row) are
// integral, possibly after scaling by a small integer, is written as
// a x <= b over integer x. Each integer column is complemented to its bound
// nearest to the LP value, x = l + x' or x = u - x' with x' >= 0, so the
// complemented LP value d_j = |x_j - bound| is >= 0. For a set S of row sides
// with sum_S b odd, (1/2) sum_S (a x <= b) rounded down is valid, and its
// violation at the LP point is
//     (1 - sum_S slack_i - sum_{j : (sum_S a_ij) odd} d_j) / 2.
// So a cut exists if a combination over GF(2) has odd right-hand side and
// total weight (slacks plus d_j of its odd columns) below 1.
//
// The mod-2 system keeps only row sides with slack < 1 and columns with d_j > 0
// (a column at its bound adds nothing to the weight); identical mod-2 rows keep
// the one with the smallest slack. Then a Gaussian elimination as in Koster,
// Zymolka and Kutschka (2009) and SCIP: columns by decreasing d_j, pivot row
// with the smallest slack among the rows containing the column, which is added
// to the other rows; after every update a row with odd right-hand side and
// weight < 1 - 2 min_violation gives a cut. Rows whose slack reaches the limit
// are dropped. Everything runs under a work budget, in a fixed order
// (deterministic).

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <numeric>

#include "milp/cuts/cut.hpp"

namespace milp {

namespace {

constexpr double kMinViolation = 0.05;  // of the zero-half cut in complemented space, before scaling
constexpr double kMaxWeight = 1.0 - 2.0 * kMinViolation;
constexpr double kMaxInt = 1e12;        // larger values are not handled (exactness of doubles)
constexpr int kMaxScale = 16;
constexpr size_t kMaxCombination = 100;  // row sides per cut
constexpr size_t kMaxCutsPerCall = 1000;   // the lightest combinations are kept
constexpr size_t kMaxCombinations = 20000;

bool integral(double v) { return std::abs(v - std::round(v)) <= 1e-9 * std::max(1.0, std::abs(v)); }

// A row side as  sum scale * sign * a_ij x_j <= rhs  (continuous columns replaced by bounds).
struct RowSide {
  Index row;
  double sign;   // +1: row <= ub, -1: -row <= -lb
  double scale;
  double rhs;    // scaled, continuous columns substituted
};

class ZeroHalfSeparator : public Separator {
 public:
  const char* name() const override { return "zerohalf"; }
  void separate(const SepaContext& ctx, std::vector<Cut>& out) override;

 private:
  void build(const MipProblem& p);
  bool built_ = false;
  SparseMatrix R_;  // row-wise problem matrix
  std::vector<RowSide> sides_;
};

void ZeroHalfSeparator::build(const MipProblem& p) {
  built_ = true;
  R_ = p.A.transpose();
  const Index m = p.num_rows();
  for (Index i = 0; i < m; ++i) {
    const Index rb = R_.begin(i), re = R_.end(i);
    if (re == rb) continue;
    bool has_int = false;
    for (Index k = rb; k < re; ++k) has_int |= p.is_int[R_.index[k]] != 0;
    if (!has_int) continue;
    const bool eq = p.row_lb[i] == p.row_ub[i];
    for (double sign : {1.0, -1.0}) {
      if (eq && sign < 0) break;  // both sides of an equation are the same mod 2
      const double side = sign > 0 ? p.row_ub[i] : -p.row_lb[i];
      if (!is_finite(side)) continue;
      // Right-hand side with continuous columns at the bound that relaxes the row.
      double rhs = side;
      bool ok = true;
      for (Index k = rb; k < re && ok; ++k) {
        const Index j = R_.index[k];
        if (p.is_int[j]) continue;
        const double a = sign * R_.value[k];
        const double bnd = a > 0 ? p.col_lb[j] : p.col_ub[j];
        ok = is_finite(bnd);
        rhs -= a * bnd;
      }
      if (!ok || std::abs(rhs) > kMaxInt) continue;
      for (int s = 1; s <= kMaxScale; ++s) {
        bool good = integral(s * rhs);
        for (Index k = rb; k < re && good; ++k) {
          const Index j = R_.index[k];
          if (!p.is_int[j]) continue;
          const double v = s * R_.value[k];
          good = integral(v) && std::abs(v) <= kMaxInt;
        }
        if (good) {
          sides_.push_back({i, sign, static_cast<double>(s), std::round(s * rhs)});
          break;
        }
      }
    }
  }
}

// A row of the mod-2 system.
struct Mod2Row {
  std::vector<Index> cols;  // sorted column ids of the system with odd coefficient
  std::vector<Index> comb;  // sorted side ids
  double slack;
  bool odd;
  bool active;
};

// Symmetric difference of two sorted vectors into `out`.
void sym_diff(const std::vector<Index>& a, const std::vector<Index>& b, std::vector<Index>& out) {
  out.clear();
  size_t i = 0, k = 0;
  while (i < a.size() && k < b.size()) {
    if (a[i] < b[k]) {
      out.push_back(a[i++]);
    } else if (b[k] < a[i]) {
      out.push_back(b[k++]);
    } else {
      ++i, ++k;
    }
  }
  out.insert(out.end(), a.begin() + static_cast<std::ptrdiff_t>(i), a.end());
  out.insert(out.end(), b.begin() + static_cast<std::ptrdiff_t>(k), b.end());
}

void ZeroHalfSeparator::separate(const SepaContext& ctx, std::vector<Cut>& out) {
  const MipProblem& p = ctx.p;
  if (!built_) build(p);
  if (sides_.empty()) return;
  const Index n = p.num_cols();
  const double* x = ctx.x;

  // Complementation: bound used per integer column (0: lower, 1: upper, -1: none) and d_j.
  std::vector<signed char> comp(n, -1);
  std::vector<double> d(n, 0.0);
  for (Index j = 0; j < n; ++j) {
    if (!p.is_int[j]) continue;
    const bool fl = is_finite(p.col_lb[j]) && std::abs(p.col_lb[j]) <= kMaxInt;
    const bool fu = is_finite(p.col_ub[j]) && std::abs(p.col_ub[j]) <= kMaxInt;
    if (fl && (!fu || x[j] - p.col_lb[j] <= p.col_ub[j] - x[j])) {
      comp[j] = 0, d[j] = std::max(0.0, x[j] - p.col_lb[j]);
    } else if (fu) {
      comp[j] = 1, d[j] = std::max(0.0, p.col_ub[j] - x[j]);
    }
  }
  constexpr double kZero = 1e-6;
  std::vector<Index> sys_col(n, -1);  // column of the mod-2 system, by decreasing d
  std::vector<Index> col_of;          // system column -> problem column
  {
    std::vector<Index> order;
    for (Index j = 0; j < n; ++j) {
      if (comp[j] >= 0 ? d[j] > kZero : p.is_int[j] != 0) order.push_back(j);
    }
    // Columns without a bound are kept with "infinite" weight (their odd sums cannot be rounded).
    std::stable_sort(order.begin(), order.end(), [&](Index a, Index b) {
      const double da = comp[a] >= 0 ? d[a] : kInf, db = comp[b] >= 0 ? d[b] : kInf;
      return da > db;
    });
    for (Index j : order) sys_col[j] = static_cast<Index>(col_of.size()), col_of.push_back(j);
  }
  auto col_weight = [&](Index c) {
    const Index j = col_of[c];
    return comp[j] >= 0 ? d[j] : kInf;
  };

  // Mod-2 rows of the row sides with slack below the limit.
  std::vector<Mod2Row> rows;
  for (Index s = 0; s < static_cast<Index>(sides_.size()); ++s) {
    const RowSide& rs = sides_[s];
    double act = 0.0, parity_rhs = rs.rhs;
    Mod2Row r;
    bool ok = true;
    for (Index k = R_.begin(rs.row); k < R_.end(rs.row) && ok; ++k) {
      const Index j = R_.index[k];
      if (!p.is_int[j]) continue;
      const double a = std::round(rs.scale * rs.sign * R_.value[k]);
      act += a * x[j];
      if (comp[j] >= 0) {
        const double bnd = comp[j] == 0 ? p.col_lb[j] : p.col_ub[j];
        const double t = a * bnd;
        ok = std::abs(t) <= kMaxInt;
        parity_rhs -= t;
      }
      if (std::fmod(std::abs(a), 2.0) == 1.0 && sys_col[j] >= 0) r.cols.push_back(sys_col[j]);
    }
    if (!ok) continue;
    r.slack = std::max(0.0, rs.rhs - act);
    if (r.slack >= kMaxWeight) continue;
    r.odd = std::fmod(std::abs(parity_rhs), 2.0) == 1.0;
    if (r.cols.empty() && !r.odd) continue;
    std::sort(r.cols.begin(), r.cols.end());
    r.comb = {s};
    r.active = true;
    rows.push_back(std::move(r));
  }
  if (rows.empty()) return;
  // Identical mod-2 rows: keep the smallest slack.
  {
    std::vector<size_t> ord(rows.size());
    std::iota(ord.begin(), ord.end(), size_t{0});
    std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {
      if (rows[a].odd != rows[b].odd) return rows[a].odd < rows[b].odd;
      if (rows[a].cols != rows[b].cols) return rows[a].cols < rows[b].cols;
      if (rows[a].slack != rows[b].slack) return rows[a].slack < rows[b].slack;
      return a < b;
    });
    std::vector<char> keep(rows.size(), 0);
    for (size_t t = 0; t < ord.size(); ++t) {
      keep[ord[t]] = t == 0 || rows[ord[t]].odd != rows[ord[t - 1]].odd || rows[ord[t]].cols != rows[ord[t - 1]].cols;
    }
    size_t w = 0;
    for (size_t k = 0; k < rows.size(); ++k) {
      if (!keep[k]) continue;
      if (w != k) rows[w] = std::move(rows[k]);
      ++w;
    }
    rows.resize(w);
  }

  std::vector<std::pair<double, std::vector<Index>>> combos;  // (weight, sides)
  auto weight = [&](const Mod2Row& r) {
    double w = r.slack;
    for (Index c : r.cols) w += col_weight(c);
    return w;
  };
  auto check = [&](const Mod2Row& r) {
    if (!r.odd) return;
    const double w = weight(r);
    if (w < kMaxWeight) combos.emplace_back(w, r.comb);
  };
  for (const Mod2Row& r : rows) check(r);

  // Gaussian elimination over the columns by decreasing d_j.
  const Index ncols = static_cast<Index>(col_of.size());
  std::vector<std::vector<Index>> col_rows(ncols);
  for (Index r = 0; r < static_cast<Index>(rows.size()); ++r) {
    for (Index c : rows[r].cols) col_rows[c].push_back(r);
  }
  const long long nnz = static_cast<long long>(p.A.index.size());
  long long work = 0;
  const long long work_limit = std::clamp<long long>(20 * nnz, 1000000LL, 50000000LL);
  std::vector<Index> tmp, cand;
  for (Index c = 0; c < ncols && work < work_limit && combos.size() < kMaxCombinations; ++c) {
    cand.clear();
    for (Index r : col_rows[c]) {
      const Mod2Row& row = rows[r];
      if (row.active && std::binary_search(row.cols.begin(), row.cols.end(), c)) cand.push_back(r);
    }
    std::sort(cand.begin(), cand.end());
    cand.erase(std::unique(cand.begin(), cand.end()), cand.end());
    col_rows[c].clear();
    work += static_cast<long long>(cand.size());
    if (cand.empty()) continue;
    Index piv = cand[0];
    for (Index r : cand) {
      const Mod2Row &a = rows[r], &b = rows[piv];
      if (a.slack < b.slack || (a.slack == b.slack && a.cols.size() < b.cols.size())) piv = r;
    }
    const Mod2Row& pr = rows[piv];
    for (Index r : cand) {
      if (r == piv) continue;
      Mod2Row& row = rows[r];
      row.slack += pr.slack;
      if (row.slack >= kMaxWeight) {
        row.active = false;
        continue;
      }
      sym_diff(row.cols, pr.cols, tmp);
      work += static_cast<long long>(row.cols.size() + pr.cols.size());
      // New columns of the row are registered in their occurrence lists.
      for (Index k : tmp) {
        if (!std::binary_search(row.cols.begin(), row.cols.end(), k)) col_rows[k].push_back(r);
      }
      row.cols.swap(tmp);
      sym_diff(row.comb, pr.comb, tmp);
      row.comb.swap(tmp);
      row.odd = row.odd != pr.odd;
      if (row.comb.size() > kMaxCombination) {
        row.active = false;
        continue;
      }
      check(row);
    }
    rows[piv].active = false;
  }
  if (combos.empty()) return;
  std::sort(combos.begin(), combos.end(), [](const auto& a, const auto& b) { return a.second < b.second; });
  combos.erase(std::unique(combos.begin(), combos.end(),
                           [](const auto& a, const auto& b) { return a.second == b.second; }),
               combos.end());
  std::stable_sort(combos.begin(), combos.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
  if (combos.size() > kMaxCutsPerCall) combos.resize(kMaxCutsPerCall);

  // Cuts in the original space.
  std::vector<double> coef(n, 0.0);
  std::vector<char> mark(n, 0);
  std::vector<Index> touched;
  for (const auto& [w, comb] : combos) {
    touched.clear();
    double rhs = 0.0;
    for (Index s : comb) {
      const RowSide& rs = sides_[s];
      rhs += rs.rhs;
      for (Index k = R_.begin(rs.row); k < R_.end(rs.row); ++k) {
        const Index j = R_.index[k];
        if (!p.is_int[j]) continue;
        if (!mark[j]) mark[j] = 1, touched.push_back(j);
        coef[j] += std::round(rs.scale * rs.sign * R_.value[k]);
      }
    }
    // Complement, check the parity, halve and round down, complement back.
    bool ok = true;
    for (Index j : touched) {
      const double a = coef[j];
      if (a == 0.0) continue;
      if (comp[j] == 0) {
        rhs -= a * p.col_lb[j];
      } else if (comp[j] == 1) {
        rhs -= a * p.col_ub[j];
      } else if (std::fmod(std::abs(a), 2.0) != 0.0) {
        ok = false;
      }
    }
    ok = ok && std::abs(rhs) <= kMaxInt && std::fmod(std::abs(rhs), 2.0) == 1.0;
    Cut cut;
    if (ok) {
      cut.origin = "zerohalf";
      double crhs = (rhs - 1.0) / 2.0;
      for (Index j : touched) {
        const double a = coef[j];
        if (a == 0.0) continue;
        if (comp[j] == 0) {
          const double q = std::floor(a / 2.0);  // q x' with x' = x - l
          if (q != 0.0) cut.idx.push_back(j), cut.val.push_back(q);
          crhs += q * p.col_lb[j];
        } else if (comp[j] == 1) {
          const double q = std::floor(-a / 2.0);  // q x' with x' = u - x
          if (q != 0.0) cut.idx.push_back(j), cut.val.push_back(-q);
          crhs -= q * p.col_ub[j];
        } else {
          cut.idx.push_back(j), cut.val.push_back(a / 2.0);
        }
      }
      cut.rhs = crhs;
    }
    for (Index j : touched) coef[j] = 0.0, mark[j] = 0;
    if (ok && !cut.idx.empty()) out.push_back(std::move(cut));
  }
}

}  // namespace

std::unique_ptr<Separator> make_zerohalf_separator() { return std::make_unique<ZeroHalfSeparator>(); }

}  // namespace milp
