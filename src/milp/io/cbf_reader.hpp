#pragma once
// Reader for the Conic Benchmark Format (CBF, versions 1-3: Friberg, "CBLIB 2014: a benchmark library for conic
// mixed-integer and continuous optimization", Math. Prog. Comp. 8, 2016), the format of CBLIB.
//
// Supported: OBJSENSE, VAR and CON with the cones F (free), L+ / L- / L= (orthants) and Q / QR (quadratic and
// rotated quadratic cones), INT, OBJACOORD / OBJBCOORD / ACOORD / BCOORD. Semidefinite (PSDVAR, PSDCON and their
// coordinates) and exponential / power cones are rejected. The model becomes:
//  - a variable cone L+ / L- / L= gives bounds; a quadratic cone (x1, ..., xk) the row x2^2 + ... + xk^2 - x1^2 <= 0
//    with x1 >= 0, a rotated one (2 x1 x2 >= x3^2 + ... + xk^2) the row x3^2 + ... + xk^2 - 2 x1 x2 <= 0 with
//    x1, x2 >= 0 (both solved as second-order cones: solve_oa's cone detection);
//  - a row cone L+ / L- / L= gives a linear row (a single-column row becomes that column's bound); a quadratic cone on
//    rows gives one column per member row (the member's own column when the row is a single column with coefficient
//    1 and no constant, else a new free column y = a x + b by an equality row) and the cone's row over these.
// Columns are named x<j>, new columns y<i>, rows r<i> (linear) and k<i> (cones). Errors are thrown as CbfError
// with "<source>:<line>: <message>". Gzip-compressed files are read transparently.

#include <stdexcept>
#include <string>

#include "milp/core/model.hpp"

namespace milp {

struct CbfError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

Model read_cbf(const std::string& path);
Model read_cbf_string(const std::string& text, const std::string& source = "<string>");

}  // namespace milp
