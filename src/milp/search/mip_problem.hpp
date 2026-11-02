#pragma once
// Internal MIP representation used by the search: always minimization,
// integer bounds already rounded.

#include <string>
#include <vector>

#include "milp/core/model.hpp"
#include "milp/core/sparse.hpp"
#include "milp/core/types.hpp"

namespace milp {

struct MipProblem {
  SparseMatrix A;  // CSC
  std::vector<double> cost, col_lb, col_ub, row_lb, row_ub;
  std::vector<char> is_int;
  double offset = 0.0;  // added to cost^T x
  double display_sense = 1.0;  // -1 if the user model maximizes (for logging only)

  Index num_cols() const { return static_cast<Index>(cost.size()); }
  Index num_rows() const { return static_cast<Index>(row_lb.size()); }

  // Build from a user model; returns false if integer bounds are inconsistent.
  static bool from_model(const Model& m, MipProblem* out);

  double objective(const std::vector<double>& x) const;
  // Max violation of bounds, rows and integrality.
  double max_violation(const std::vector<double>& x, double* int_viol = nullptr) const;
};

}  // namespace milp
