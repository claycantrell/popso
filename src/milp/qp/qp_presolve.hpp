#pragma once
// Presolve for problems with a quadratic objective (QP and MIQP)
//
//   min c^T x + 1/2 x^T Q x + offset   s.t.  row_lb <= A x <= row_ub,  lb <= x <= ub
//
// The MILP presolve (presolve/presolve.hpp) assumes a linear objective; this one applies only
// reductions that stay valid with Q and keep the dual information recoverable:
//   - fixed columns: removed; their Hessian entries move into the costs of the other columns
//     (c_k += Q_kj v) and into the offset (c_j v + 1/2 Q_jj v^2), their rows' bounds are shifted;
//   - empty columns (no row entries and no off-diagonal Hessian entries): fixed at the minimizer
//     of c_j x + 1/2 Q_jj x^2 over their bounds (the better rounding for an integer column);
//   - columns without Hessian entries: dual fixing by locks (as for an LP column);
//   - empty rows: removed (or the problem is infeasible);
//   - singleton rows: turned into column bounds (rounded for integer columns);
//   - activity: rows that cannot be violated within the column bounds are removed (one side
//     at a time), rows that cannot be satisfied prove infeasibility; for MIQPs the bounds of
//     integer columns are tightened from the row activities.
// No column is substituted or aggregated (Q would have to be transformed), so the postsolve of
// x only restores the removed columns' values. Duals (continuous problems): the reduced
// problem's row duals are kept, removed rows get 0 except singleton rows whose bound is active
// in the final bounds of their column (the column's reduced cost moves to the row), and all
// reduced costs are recomputed from the original data, d = c + Q x - A^T y.

#include <vector>

#include "milp/core/log.hpp"
#include "milp/core/types.hpp"
#include "milp/search/mip_problem.hpp"

namespace milp {

enum class QpPresolveStatus { Reduced, Infeasible };

struct QpPresolveStats {
  Index rows_before = 0, cols_before = 0, rows_after = 0, cols_after = 0;
  long long fixed_cols = 0, empty_cols = 0, dual_fixed = 0;
  long long empty_rows = 0, singleton_rows = 0, redundant_rows = 0, bound_changes = 0;
  int rounds = 0;
};

class QpPresolver {
 public:
  // p: minimization form (MipProblem::from_model); p.implications must be empty.
  QpPresolveStatus run(const MipProblem& p, const Logger& log);

  const MipProblem& reduced() const { return reduced_; }
  const QpPresolveStats& stats() const { return stats_; }
  // Original column of each reduced column / original row of each reduced row.
  const std::vector<Index>& col_map() const { return col_map_; }
  const std::vector<Index>& row_map() const { return row_map_; }

  // x of the reduced problem -> x of the original problem.
  std::vector<double> postsolve(const std::vector<double>& x_reduced) const;
  // Original x and the reduced problem's row duals y_reduced (min form: d = c + Q x - A^T y)
  // -> row duals y and reduced costs d of the original problem.
  void postsolve_duals(const std::vector<double>& x, const std::vector<double>& y_reduced, std::vector<double>& y,
                       std::vector<double>& d) const;
  // Original-space vector (e.g. a start, NaN allowed) restricted to the reduced columns.
  std::vector<double> crush(const std::vector<double>& x) const;

 private:
  MipProblem orig_;  // copy of the original problem (dual postsolve)
  MipProblem reduced_;
  QpPresolveStats stats_;
  std::vector<Index> col_map_, row_map_;
  std::vector<double> fixed_val_;  // per original column: value if removed
  std::vector<char> col_removed_;
  // Singleton rows in removal order: row, column, coefficient.
  struct Singleton {
    Index row, col;
    double a;
  };
  std::vector<Singleton> singletons_;
  // Per original column: the singleton row whose bound is the column's final lower / upper bound (-1: none).
  std::vector<Index> lb_row_, ub_row_;
};

}  // namespace milp
