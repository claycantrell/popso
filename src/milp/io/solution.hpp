#pragma once
// Solution files in the common "name value" .sol format:
//
//   # Solution for model <name>
//   # Status = OPTIMAL
//   # Objective value = <17 significant digits>
//   <varname> <value>
//   ...
//
// The objective line and the variable values are only written if the result
// has a solution. Values use the shortest exact round-trip representation.

#include <string>
#include <vector>

#include "milp/core/model.hpp"
#include "milp/solver/solver.hpp"

namespace milp {

std::string format_solution(const Model& model, const SolveResult& result);

// Throws std::runtime_error on I/O failure.
void write_solution(const std::string& path, const Model& model, const SolveResult& result);

// Reads "<name> <value>" lines (lines starting with '#' and blank lines are
// skipped) and maps them onto the model's columns. Columns not mentioned in
// the file are NaN (partial solutions are allowed, e.g. for warm starts).
// Unknown names, malformed lines or bad numbers throw std::runtime_error.
std::vector<double> read_solution(const std::string& path, const Model& model);
std::vector<double> read_solution_string(const std::string& text, const Model& model);

}  // namespace milp
