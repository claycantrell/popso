#include "milp/solver/postsolve.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

#include "milp/lp/simplex.hpp"

namespace milp {

bool polish_point(const MipProblem& p, const Timer& timer, double time_limit, std::vector<double>& x,
                  double primal_tol) {
  std::vector<double> lb = p.col_lb, ub = p.col_ub;
  bool any_cont = false;
  for (Index j = 0; j < p.num_cols(); ++j) {
    if (p.is_int[j]) lb[j] = ub[j] = std::round(x[j]);
    else any_cont = true;
  }
  std::vector<double> y = x;
  if (any_cont) {
    DualSimplex lp;
    lp.load(p.A, p.cost, lb, ub, p.row_lb, p.row_ub);
    LpOptions opt;
    opt.perturb = false;
    opt.tol.primal_feas = primal_tol;
    opt.time_limit = std::max(1.0, time_limit - timer.elapsed());
    if (lp.solve(opt) != LpStatus::Optimal) return false;
    y.assign(lp.x(), lp.x() + p.num_cols());
  } else {
    for (Index j = 0; j < p.num_cols(); ++j) y[j] = lb[j];
  }
  if (!(p.max_violation(y) <= 1e-6)) return false;
  x = std::move(y);
  return true;
}

bool postsolve_result(const MipProblem& p, const Presolver& pre, const SolverParams& params, const Timer& timer,
                      const Logger& log, MipResult* mr) {
  if (!mr->has_solution) return true;
  mr->x = pre.postsolve(mr->x);
  const double viol = p.max_violation(mr->x);
  const double scale = std::max(1.0, std::abs(mr->objective));
  if (!(viol <= 1e-6)) {
    // The polished point rounds the integers exactly, so its objective may differ from the
    // presolved one by round-off (n5-3: 8429 vs 8428.999991). Accept it within a relative 1e-6,
    // and at any objective when the search stopped at a limit (a feasible point beats none).
    std::vector<double> fixed = mr->x;
    if (!polish_point(p, timer, params.time_limit, fixed) ||
        (mr->status == Status::Optimal && !(p.objective(fixed) <= mr->objective + 1e-6 * scale))) {
      log.log(1, "WARNING: postsolved solution violates the original problem by %.3g; re-solving without presolve\n",
              viol);
      mr->has_solution = false;
      mr->x.clear();
      return false;
    }
    log.log(1, "Postsolved solution violated the original problem by %.3g; repaired by an LP over the continuous "
            "columns\n", viol);
    mr->x = std::move(fixed);
  } else if (const double obj = p.objective(mr->x); obj > mr->objective + 1e-9 * scale) {
    // Postsolve keeps the objective up to round-off, except where it clamps a substituted column into
    // its bounds. The reduced point may violate the row that holds such a bound within the feasibility
    // tolerance, and the clamp moves that violation from the bound into the column's objective term:
    // with a large cost a lot (binkar10_1: penalty columns with cost 1e6 substituted out of their
    // equations, 2.5e-8 below their lower bound 0 in an LP point; the clamped point came out 0.0256
    // worse than the search's incumbent, beyond the gap 1e-6, and was reported as optimal). Re-optimize
    // the continuous columns, at a tight tolerance: that LP has the same costs.
    for (const double tol : {1e-9, 1e-7}) {
      std::vector<double> fixed = mr->x;
      if (!polish_point(p, timer, params.time_limit, fixed, tol)) continue;
      if (p.objective(fixed) < obj) {
        log.log(2, "Postsolved solution: objective %.10g against %.10g in the presolved problem; continuous columns "
                "re-optimized: %.10g\n", p.display_sense * obj, p.display_sense * mr->objective,
                p.display_sense * p.objective(fixed));
        mr->x = std::move(fixed);
      }
      break;
    }
  }
  const double obj = p.objective(mr->x);
  if (obj > mr->objective + 1e-6 * scale) {
    // Still worse: the search's incumbent itself was better than its integer assignment allows, by
    // what the feasibility tolerance buys at large costs (binkar10_1 seeds 1 and 5 deterministic:
    // 6742.172882 in the presolved problem, re-optimized 6742.200024, the optimum). Its cutoff was
    // that much too optimistic. Re-solving without presolve to certify the gap cost 3-4x the time
    // there, for the same answer; not done.
    log.log(2, "Postsolved objective %.10g is %.3g worse than the search's %.10g\n", p.display_sense * obj,
            obj - mr->objective, p.display_sense * mr->objective);
  }
  mr->objective = obj;
  mr->bound = std::min(mr->bound, obj);
  return true;
}

}  // namespace milp
