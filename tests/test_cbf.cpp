// The CBF reader (io/cbf_reader.hpp): cones on variables and on rows, integers, bounds from singleton rows, errors.

#include <cmath>
#include <string>

#include "milp/io/cbf_reader.hpp"
#include "milp/solver/solver.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

SolveResult solve_cbf(const std::string& text) {
  SolverParams p;
  p.verbosity = 0;
  p.threads = 1;
  p.mip_rel_gap = 1e-9;
  p.mip_abs_gap = 1e-9;
  return solve(read_cbf_string(text), p);
}

}  // namespace

TEST(cbf_cones_on_variables) {
  // min t s.t. (t, x, y) in Q^3, x + y = 2 (an L= row): sqrt(2).
  const std::string text =
      "# comment\nVER\n3\n\nOBJSENSE\nMIN\n\nVAR\n3 1\nQ 3\n\nCON\n1 1\nL= 1\n\n"
      "OBJACOORD\n1\n0 1.0\n\nACOORD\n2\n0 1 1.0\n0 2 1.0\n\nBCOORD\n1\n0 -2.0\n";
  const Model m = read_cbf_string(text);
  CHECK(m.num_cols() == 3);
  CHECK(m.num_rows() == 2);  // the L= row and the cone row
  CHECK(m.col_lb[0] == 0.0);  // the cone's head
  const SolveResult r = solve_cbf(text);
  CHECK(r.status == Status::Optimal);
  CHECK_NEAR(r.objective, std::sqrt(2.0), 1e-7);
}

TEST(cbf_rotated_cone_on_rows_and_integers) {
  // max -(u + v) s.t. 2 u v >= (x - 1)^2 written on rows (u, v, x - 1) in QR^3, x >= 3.5 by a singleton L+ row,
  // x integer: x = 4, u v >= 4.5, min u + v = 2 sqrt(4.5) = 4.2426...
  const std::string text =
      "VER\n1\nOBJSENSE\nMAX\nVAR\n3 1\nF 3\nINT\n1\n2\nCON\n4 2\nQR 3\nL+ 1\n"
      "OBJACOORD\n2\n0 -1\n1 -1\nACOORD\n4\n0 0 1\n1 1 1\n2 2 1\n3 2 1\nBCOORD\n2\n2 -1\n3 -3.5\n";
  const Model m = read_cbf_string(text);
  CHECK(m.col_lb[2] == 4.0);  // x >= 3.5 as a bound, rounded up for an integer
  CHECK(m.col_type[2] == VarType::Integer);
  CHECK(m.sense == ObjSense::Maximize);
  const SolveResult r = solve_cbf(text);
  CHECK(r.status == Status::Optimal);
  CHECK_NEAR(r.objective, -2.0 * std::sqrt(4.5), 1e-6);
}

TEST(cbf_unsupported_and_invalid) {
  bool caught = false;
  try {
    read_cbf_string("VER\n3\nVAR\n2 1\nEXP 2\n");
  } catch (const CbfError& e) {
    caught = std::string(e.what()).find("EXP") != std::string::npos;
  }
  CHECK(caught);
  caught = false;
  try {
    read_cbf_string("VER\n3\nPSDVAR\n1\n2\n");
  } catch (const CbfError& e) {
    caught = std::string(e.what()).find("PSDVAR") != std::string::npos;
  }
  CHECK(caught);
  caught = false;
  try {
    read_cbf_string("VER\n3\nVAR\n2 1\nF 3\n");  // the cones cover 3 of 2
  } catch (const CbfError&) {
    caught = true;
  }
  CHECK(caught);
}
