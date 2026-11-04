// Gomory mixed-integer cuts from simplex tableau rows.
//
// The tableau row of a basic integer variable x_B reads
//   x_B + sum_{j nonbasic} alpha_j x_j = 0      (over structurals and logicals).
// Nonbasic variables are shifted to their active bound (complemented at upper),
// y_j = x_j - l_j or y_j = u_j - x_j >= 0, giving x_B + sum a_j y_j = beta, from
// which the GMI cut sum g_j y_j >= 1 follows. Logicals s_i = -a_i x are then
// substituted back to obtain a cut over structural variables only.

#include <algorithm>
#include <cmath>

#include "milp/cuts/cut.hpp"

namespace milp {

namespace {

class GomorySeparator : public Separator {
 public:
  const char* name() const override { return "gomory"; }
  void separate(const SepaContext& ctx, std::vector<Cut>& out) override;

 private:
  std::vector<double> alpha_, coef_;
  std::vector<Index> anz_;  // nonzero positions of alpha_ (ascending)
  std::vector<char> int_row_, mark_;
};

bool near_int(double v, double tol) { return std::abs(v - std::round(v)) <= tol; }

// Integer, or implied integer (MipProblem::implied_int: integral in some optimal solution) with integral
// global bounds.
bool int_col(const MipProblem& p, Index j) {
  if (p.is_int[j]) return true;
  auto ib = [](double b) { return !is_finite(b) || b == std::round(b); };
  return p.is_implied_int(j) && ib(p.col_lb[j]) && ib(p.col_ub[j]);
}

void GomorySeparator::separate(const SepaContext& ctx, std::vector<Cut>& out) {
  const MipProblem& p = ctx.p;
  DualSimplex& lp = ctx.lp;
  const Index n = lp.num_cols(), m = lp.num_rows();
  const SparseMatrix& R = lp.row_matrix();
  const double* x = ctx.x;

  // Rows whose logical is integer for integer x: integer coefficients on integer columns.
  int_row_.assign(m, 1);
  for (Index i = 0; i < m; ++i) {
    for (Index k = R.begin(i); k < R.end(i) && int_row_[i]; ++k) {
      int_row_[i] = int_col(p, R.index[k]) && R.value[k] == std::round(R.value[k]);
    }
  }

  // Candidate rows: fractional basic integer structurals, most fractional first.
  std::vector<std::pair<double, Index>> cand;
  for (Index j : lp.basic_vars()) {
    if (j >= n || !p.is_int[j]) continue;
    const double f = x[j] - std::floor(x[j]);
    if (f < 0.01 || f > 0.99) continue;
    cand.emplace_back(std::abs(f - 0.5), j);
  }
  std::sort(cand.begin(), cand.end());
  const size_t max_rows = 500;
  if (cand.size() > max_rows) cand.resize(max_rows);

  coef_.assign(n, 0.0);
  mark_.assign(n, 0);
  alpha_.assign(n + m, 0.0);
  anz_.clear();
  // Rows with more than twice the pool's density limit in nonbasic structurals give cuts the pool
  // rejects (cancellation through the logicals is rare): skip them before the second pass
  // (neos-960392: ~40% of the rows have 50k+ of 58k columns).
  const Index max_struct = ctx.max_nnz > 0 ? 2 * ctx.max_nnz : n + 1;
  std::vector<Index> nz;
  for (const auto& [score, jb] : cand) {
    for (Index j : anz_) alpha_[j] = 0.0;
    anz_.clear();
    if (!lp.tableau_row(jb, alpha_, anz_)) continue;
    // First pass: beta and validity of the row (positions outside anz_ are zero).
    double beta = 0.0, amax = 0.0;
    bool ok = true;
    Index nstruct = 0;
    for (Index j : anz_) {
      if (lp.status(j) == VarStatus::Basic) continue;
      const double a = alpha_[j];
      if (std::abs(a) < 1e-11) continue;  // numerically zero tableau entry
      amax = std::max(amax, std::abs(a));
      const VarStatus st = lp.status(j);
      const double l = lp.var_lower(j), u = lp.var_upper(j);
      if (l == u) {
        beta -= a * l;
        continue;
      }
      if (j < n && ++nstruct > max_struct) {
        ok = false;  // too dense
        break;
      }
      if (st == VarStatus::AtZero) {
        ok = false;  // free nonbasic variable
        break;
      }
      const double b = st == VarStatus::AtUpper ? u : l;
      if (!is_finite(b)) ok = false;
      beta -= a * b;
    }
    if (!ok || amax > 1e6) continue;
    if (std::abs(beta - x[jb]) > 1e-6 * (1.0 + std::abs(x[jb]))) continue;
    const double f0 = beta - std::floor(beta);
    if (f0 < 0.01 || f0 > 0.99) continue;

    // Second pass: GMI coefficients g_j for y_j, accumulated in x space.
    //   sum_j g_j y_j >= 1,  y_j = sign_j (x_j - bound_j).
    double rhs = 1.0;  // of the >= cut in x space
    auto add = [&](Index j, double v) {
      if (!mark_[j]) mark_[j] = 1, nz.push_back(j);
      coef_[j] += v;
    };
    for (Index j : anz_) {
      if (lp.status(j) == VarStatus::Basic) continue;
      const double a0 = alpha_[j];
      if (std::abs(a0) < 1e-11) continue;
      const double l = lp.var_lower(j), u = lp.var_upper(j);
      if (l == u) continue;
      const bool at_upper = lp.status(j) == VarStatus::AtUpper;
      const double sign = at_upper ? -1.0 : 1.0;
      const double bound = at_upper ? u : l;
      const double a = sign * a0;  // coefficient of y_j
      bool integer;
      if (j < n) {
        integer = p.is_int[j] || (p.is_implied_int(j) && bound == std::round(bound));
      } else {
        integer = int_row_[j - n] && near_int(bound, 1e-9);
      }
      double g;
      if (integer) {
        const double fj = a - std::floor(a);
        g = fj <= f0 ? fj / f0 : (1.0 - fj) / (1.0 - f0);
      } else {
        g = a > 0 ? a / f0 : -a / (1.0 - f0);
      }
      if (g == 0.0) continue;
      // g * y_j = g * sign * x_j - g * sign * bound.
      rhs += g * sign * bound;
      const double gx = g * sign;
      if (j < n) {
        add(j, gx);
      } else {
        // x_j = s_i = -a_i x.
        const Index i = j - n;
        for (Index k = R.begin(i); k < R.end(i); ++k) add(R.index[k], -gx * R.value[k]);
      }
    }
    Cut c;
    c.origin = "gomory";
    for (Index k : nz) {
      c.idx.push_back(k);
      c.val.push_back(-coef_[k]);
      coef_[k] = 0.0, mark_[k] = 0;
    }
    nz.clear();
    c.rhs = -rhs;
    out.push_back(std::move(c));
  }
}

}  // namespace

std::unique_ptr<Separator> make_gomory_separator() { return std::make_unique<GomorySeparator>(); }

}  // namespace milp
