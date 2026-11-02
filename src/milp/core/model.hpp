#pragma once
// The user-facing MILP model:
//
//   min/max  c^T x + offset
//   s.t.     row_lb <= A x <= row_ub
//            col_lb <=   x <= col_ub
//            x_j integer for j in integer columns
//
// Rows are ranged; equality rows have row_lb == row_ub, free rows (other
// than the objective) are allowed. The constraint matrix is stored by column
// (CSC). Use ModelBuilder to assemble a model incrementally.

#include <string>
#include <vector>

#include "milp/core/sparse.hpp"
#include "milp/core/types.hpp"

namespace milp {

struct Model {
  std::string name;
  ObjSense sense = ObjSense::Minimize;
  double obj_offset = 0.0;

  // Columns.
  std::vector<std::string> col_names;
  std::vector<double> obj;
  std::vector<double> col_lb;
  std::vector<double> col_ub;
  std::vector<VarType> col_type;

  // Rows.
  std::vector<std::string> row_names;
  std::vector<double> row_lb;
  std::vector<double> row_ub;

  SparseMatrix A;  // CSC, num_major = num_cols, num_minor = num_rows

  Index num_cols() const { return static_cast<Index>(obj.size()); }
  Index num_rows() const { return static_cast<Index>(row_lb.size()); }
  Index num_nonzeros() const { return A.nnz(); }
  Index num_integers() const;
  bool is_mip() const { return num_integers() > 0; }

  // Objective value of x in the original sense (including offset).
  double objective_value(const std::vector<double>& x) const;

  // Maximum absolute violation of bounds, rows and (optionally) integrality.
  struct Violation {
    double bound = 0.0, row = 0.0, integrality = 0.0;
    double max() const;
  };
  Violation violation(const std::vector<double>& x) const;

  // Throws std::runtime_error on inconsistent dimensions / NaNs.
  void validate() const;
};

// Incremental builder; rows and columns may be added in any order, nonzeros
// are collected as triplets and assembled by build().
class ModelBuilder {
 public:
  Index add_col(std::string name, double obj, double lb, double ub, VarType type = VarType::Continuous);
  Index add_row(std::string name, double lb, double ub);
  void add_coef(Index row, Index col, double value);
  Model& model() { return m_; }
  Model build();

 private:
  Model m_;
  std::vector<Index> trip_col_, trip_row_;
  std::vector<double> trip_val_;
};

}  // namespace milp
