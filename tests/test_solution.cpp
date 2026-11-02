// Tests for the solution writer / reader.

#include <cmath>
#include <filesystem>
#include <stdexcept>

#include "milp/io/solution.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

Model three_vars() {
  ModelBuilder b;
  b.model().name = "sol_model";
  b.add_col("x", 1.0, 0.0, 10.0);
  b.add_col("y_long_name[1,2]", 1.0, -kInf, kInf);
  b.add_col("z", 1.0, 0.0, 1.0, VarType::Binary);
  return b.build();
}

}  // namespace

TEST(solution_format_optimal) {
  Model m = three_vars();
  SolveResult r;
  r.status = Status::Optimal;
  r.has_solution = true;
  r.objective = 0.1 + 0.2;
  r.x = {1.0, -2.5, 0.0};
  const std::string text = format_solution(m, r);
  CHECK_EQ(text,
           "# Solution for model sol_model\n"
           "# Status = OPTIMAL\n"
           "# Objective value = 0.30000000000000004\n"
           "x 1\n"
           "y_long_name[1,2] -2.5\n"
           "z 0\n");
}

TEST(solution_format_no_solution) {
  Model m = three_vars();
  SolveResult r;
  r.status = Status::Infeasible;
  CHECK_EQ(format_solution(m, r), "# Solution for model sol_model\n# Status = INFEASIBLE\n");
  r.status = Status::TimeLimit;
  CHECK_EQ(format_solution(m, r), "# Solution for model sol_model\n# Status = TIME_LIMIT\n");
}

TEST(solution_roundtrip_exact) {
  Model m = three_vars();
  SolveResult r;
  r.status = Status::TimeLimit;
  r.has_solution = true;
  r.x = {1.0 / 3.0, -1e-300, 0.1 + 0.2};
  r.objective = m.objective_value(r.x);
  const auto path = std::filesystem::temp_directory_path() / "milp_test_sol.sol";
  write_solution(path.string(), m, r);
  std::vector<double> x = read_solution(path.string(), m);
  REQUIRE(x.size() == 3);
  CHECK(x == r.x);  // bit-exact
  std::filesystem::remove(path);
}

TEST(solution_read_partial_and_errors) {
  Model m = three_vars();
  std::vector<double> x = read_solution_string("# comment\n\n  z   1 \r\nx 2.5\n", m);
  REQUIRE(x.size() == 3);
  CHECK_EQ(x[0], 2.5);
  CHECK(std::isnan(x[1]));
  CHECK_EQ(x[2], 1.0);
  CHECK_THROWS(read_solution_string("w 1\n", m), std::runtime_error);
  CHECK_THROWS(read_solution_string("x one\n", m), std::runtime_error);
  CHECK_THROWS(read_solution_string("x 1 2\n", m), std::runtime_error);
  CHECK_THROWS(read_solution("/nonexistent/dir/x.sol", m), std::runtime_error);
  CHECK_THROWS(write_solution("/nonexistent/dir/x.sol", m, SolveResult{}), std::runtime_error);
}
