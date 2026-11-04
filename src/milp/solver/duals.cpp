#include "milp/solver/duals.hpp"

#include <algorithm>
#include <cmath>

#include "milp/core/interrupt.hpp"
#include "milp/core/timer.hpp"
#include "milp/lp/crossover.hpp"

namespace milp {

LpDualResult recover_lp_duals(const MipProblem& p, const std::vector<double>& col_lb_in,
                              const std::vector<double>& col_ub_in, const std::vector<double>& x, double time_limit) {
  LpDualResult r;
  const Index n = p.num_cols(), m = p.num_rows();
  const std::vector<double>& lb = col_lb_in.empty() ? p.col_lb : col_lb_in;
  const std::vector<double>& ub = col_ub_in.empty() ? p.col_ub : col_ub_in;
  if (static_cast<Index>(x.size()) != n) return r;
  if (m == 0) {
    // No rows: d = c, each column at the bound its cost points to.
    r.row_dual.clear();
    r.reduced_cost = p.cost;
    r.col_status.resize(static_cast<size_t>(n));
    r.objective = p.offset;
    for (Index j = 0; j < n; ++j) {
      const double c = p.cost[j];
      VarStatus st = VarStatus::AtZero;
      if (lb[j] == ub[j]) st = VarStatus::AtLower;
      else if (c > 0.0) st = VarStatus::AtLower;
      else if (c < 0.0) st = VarStatus::AtUpper;
      else if (is_finite(lb[j])) st = VarStatus::AtLower;
      else if (is_finite(ub[j])) st = VarStatus::AtUpper;
      if ((st == VarStatus::AtLower && !is_finite(lb[j])) || (st == VarStatus::AtUpper && !is_finite(ub[j]))) return r;
      r.col_status[j] = st;
      r.objective += c * (st == VarStatus::AtLower ? lb[j] : st == VarStatus::AtUpper ? ub[j] : 0.0);
    }
    r.ok = true;
    return r;
  }
  Timer timer;
  // The crossover honours interrupts and this deadline; an interrupted recovery gives no duals.
  DeadlineScope deadline(time_limit);
  if (interrupted_or_expired()) return r;
  DualSimplex lp;
  lp.load(p.A, p.cost, lb, ub, p.row_lb, p.row_ub);
  std::vector<double> xc(x);
  for (Index j = 0; j < n; ++j) xc[j] = std::clamp(xc[j], lb[j], ub[j]);
  const Basis b = crossover(p.A, p.cost, lb, ub, p.row_lb, p.row_ub, xc, {});
  if (interrupted_or_expired()) return r;
  lp.set_basis(b);
  LpOptions lo;
  lo.primal_first = true;
  lo.time_limit = std::max(0.0, time_limit - timer.elapsed());
  LpStatus st = lp.solve(lo);
  if (st != LpStatus::Optimal && st != LpStatus::TimeLimit && !interrupted()) {
    // Numerical trouble from the crossover basis: solve from the slack basis instead.
    lp.reset_basis();
    lo.primal_first = false;
    lo.time_limit = std::max(0.0, time_limit - timer.elapsed());
    st = lp.solve(lo);
  }
  r.iterations = lp.iterations();
  if (st != LpStatus::Optimal) return r;
  r.row_dual = lp.row_duals();
  r.reduced_cost = lp.reduced_costs();
  r.objective = lp.objective() + p.offset;
  r.col_status.resize(static_cast<size_t>(n));
  r.row_status.resize(static_cast<size_t>(m));
  for (Index j = 0; j < n; ++j) r.col_status[j] = lp.status(j);
  for (Index i = 0; i < m; ++i) {
    // The logical s_i = -a_i x: at its lower bound -row_ub the row is at its upper bound.
    const VarStatus s = lp.status(n + i);
    if (s == VarStatus::Basic || s == VarStatus::AtZero) r.row_status[i] = s;
    else if (p.row_lb[i] == p.row_ub[i]) r.row_status[i] = VarStatus::AtLower;  // equality row
    else r.row_status[i] = s == VarStatus::AtLower ? VarStatus::AtUpper : VarStatus::AtLower;
  }
  r.ok = true;
  return r;
}

DualObjectiveCheck check_dual_objective(double basis_objective, double solution_objective, bool fixed_mip) {
  const double tol = 1e-6 * (1.0 + std::abs(solution_objective));
  const double diff = basis_objective - solution_objective;
  if (std::abs(diff) <= tol) return DualObjectiveCheck::Consistent;
  if (fixed_mip && diff < -tol) return DualObjectiveCheck::BetterThanSolution;
  return DualObjectiveCheck::Mismatch;  // also NaN
}

}  // namespace milp
