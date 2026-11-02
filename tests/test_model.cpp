// Tests for Model and ModelBuilder.

#include <stdexcept>

#include "milp/core/model.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

// min x + 2y - z + 3   s.t.  x + y >= 1,  y - z <= 2,  0 <= x <= 4, y int, z binary
Model small_model() {
  ModelBuilder b;
  Index x = b.add_col("x", 1.0, 0.0, 4.0);
  Index y = b.add_col("y", 2.0, -kInf, kInf, VarType::Integer);
  Index z = b.add_col("z", -1.0, -5.0, 7.0, VarType::Binary);
  Index r0 = b.add_row("r0", 1.0, kInf);
  Index r1 = b.add_row("r1", -kInf, 2.0);
  b.add_coef(r0, x, 1.0);
  b.add_coef(r0, y, 1.0);
  b.add_coef(r1, y, 1.0);
  b.add_coef(r1, z, -1.0);
  b.model().obj_offset = 3.0;
  return b.build();
}

}  // namespace

TEST(model_builder_dimensions) {
  Model m = small_model();
  CHECK_EQ(m.num_cols(), 3);
  CHECK_EQ(m.num_rows(), 2);
  CHECK_EQ(m.num_nonzeros(), 4);
  CHECK_EQ(m.num_integers(), 2);
  CHECK(m.is_mip());
  CHECK_EQ(m.col_names[1], "y");
  CHECK_EQ(m.row_names[1], "r1");
  // Binary bounds are clamped to [0, 1].
  CHECK_EQ(m.col_lb[2], 0.0);
  CHECK_EQ(m.col_ub[2], 1.0);
  // Column y holds rows 0 and 1.
  CHECK_EQ(m.A.end(1) - m.A.begin(1), 2);
}

TEST(model_builder_resets) {
  ModelBuilder b;
  b.add_col("a", 1.0, 0.0, 1.0);
  Model m1 = b.build();
  Model m2 = b.build();
  CHECK_EQ(m1.num_cols(), 1);
  CHECK_EQ(m2.num_cols(), 0);
  CHECK(!m2.is_mip());
}

TEST(model_objective_value) {
  Model m = small_model();
  CHECK_NEAR(m.objective_value({1.0, 2.0, 1.0}), 1.0 + 4.0 - 1.0 + 3.0, 1e-12);
}

TEST(model_violation) {
  Model m = small_model();
  auto v = m.violation({1.0, 0.0, 1.0});
  CHECK_EQ(v.max(), 0.0);
  // x below row 0, y fractional, x above bound.
  v = m.violation({5.0, 0.5, 0.0});
  CHECK_NEAR(v.bound, 1.0, 1e-12);
  CHECK_NEAR(v.integrality, 0.5, 1e-12);
  CHECK_NEAR(v.row, 0.0, 1e-12);
  v = m.violation({0.0, 0.0, 0.0});
  CHECK_NEAR(v.row, 1.0, 1e-12);
  v = m.violation({0.0, 4.0, 0.0});
  CHECK_NEAR(v.row, 2.0, 1e-12);
  CHECK_NEAR(v.max(), 2.0, 1e-12);
}

TEST(model_validate_errors) {
  Model m = small_model();
  m.col_lb.pop_back();
  CHECK_THROWS(m.validate(), std::runtime_error);
  Model m2 = small_model();
  m2.obj[0] = std::nan("");
  CHECK_THROWS(m2.validate(), std::runtime_error);
  Model m3 = small_model();
  m3.row_names.push_back("extra");
  CHECK_THROWS(m3.validate(), std::runtime_error);
}
