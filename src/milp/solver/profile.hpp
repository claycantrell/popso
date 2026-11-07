#pragma once
// Model and run profiles (popso --profile): JSON with the structure of a model, its numerics and
// how a solve went, for users who cannot share a model at all. Everything is aggregated: no
// names, no per-row or per-column data, no file paths. From profiles we learn which structures
// users solve and how the solver behaves on them (heuristic effort per class, look-alike test
// models).
//
// Schema ("profile_version": 1; null stands for a missing or infinite number):
//   profile_version
//   solver:   name, version, threads (used), time_limit, deterministic, seed, racing, presolve,
//             mip_rel_gap, lp_method
//   model:    sense ("minimize" / "maximize"), rows, cols, nonzeros, objective_nonzeros,
//             col_types {continuous, binary, integer, semi_continuous, semi_integer},
//             col_bounds {free, lower, upper, boxed, fixed},
//             row_senses {equal, less, greater, ranged, free},
//             row_classes {empty, free, singleton, aggregation, precedence, variable_bound,
//                          set_partitioning, set_packing, set_covering, cardinality,
//                          invariant_knapsack, equation_knapsack, bin_packing, knapsack,
//                          integer_knapsack, mixed_binary, general_linear}  (see RowClass),
//             rows_continuous_only (rows whose columns are all continuous or semi-continuous;
//             they count as mixed_binary too),
//             sos1, sos2, sos_members, indicators, indicator_nonzeros,
//             quadratic {nonzeros, diagonal, cols}  (Q's lower triangle; cols: columns in Q)
//   numerics: matrix, objective, bounds, rhs, quadratic: {min, max} of the nonzero absolute
//             values (bounds and rhs: finite ones; the indicator rows are not included),
//             rows_coef_range_above_1e6 (max |a| / min |a| of the row above 1e6),
//             rows_with_tiny_entries (an entry below 1e-9 of the row's largest |a|),
//             row_length, col_length: {min, median, max} (entries per row / column)
//   presolve: ran, rows_before, cols_before, nonzeros_before, rows, cols, nonzeros, time
//             (SolveStats::presolve: the problem after the logical reformulation, first presolve)
//   run:      status, has_solution, objective, bound, gap (|objective - bound| / max(|objective|,
//             1e-10), as popso prints it), time, nodes, lp_iterations, barrier_iterations,
//             root_bound, restarts, tree_restarts,
//             incumbents: [{time, objective, source, count, last_time, last_objective}]
//             (SolveResult::incumbent_trace: a run of quick improvements from one source is one
//             entry; empty unless SolverParams::trace_incumbents was set)

#include <cstdint>
#include <string>
#include <vector>

#include "milp/core/model.hpp"
#include "milp/solver/solver.hpp"

namespace milp {

// Row classes of the MIPLIB 2017 constraint classification (Gleixner et al., "MIPLIB 2017",
// Math. Prog. Comp. 2021, the table of constraint types), in the table's order; classify_rows
// gives each row the first class that matches. Exact rules:
//  - Column kinds: binary = an integer column (Integer or Binary type, not semi-integer) whose
//    bounds admit only 0 and 1 (lb > -1, ub < 2); integer = the other integer and semi-integer
//    columns; continuous = continuous and semi-continuous columns.
//  - The row's entries are its nonzeros (n of them); its sides lb <= a^T x <= ub. An equation has
//    lb == ub. An inequality is used in "<= form": a^T x <= ub if ub is finite, -a^T x <= -lb if lb
//    is finite (a ranged row has both forms; a class matches if it matches one of them). In a
//    <= form c^T x <= beta over binaries, a column with c_k < 0 is complemented (x = 1 - y), which
//    gives beta' = beta + sum of |c_k| over c_k < 0 and all coefficients |c_k|; equations likewise
//    (b' = b + sum of |a_k| over a_k < 0). "Equal" and "integral" allow 1e-9 (relative above 1).
//  1 empty: n = 0.                       2 free: no finite side.
//  3 singleton: n = 1.                   4 aggregation: n = 2, equation.
//  5 precedence: n = 2, inequality, a_1 = -a_2, both columns of the same kind.
//  6 variable_bound: n = 2, inequality, at least one binary column.
//  Classes 7-11 need all columns binary and all |a_k| equal (to s):
//  7 set_partitioning: equation, b' / s = 1.
//  8 set_packing: a <= form with beta' / s = 1.
//  9 set_covering: inequality with finite lb, as ">= form" sum of literals >= 1:
//    (lb + s * #{a_k < 0}) / s = 1.
// 10 cardinality: equation, b' / s integral and >= 2.
// 11 invariant_knapsack: a <= form with beta' / s integral and >= 2.
//  Classes 12-14 need all columns binary and all a_k integral:
// 12 equation_knapsack: equation, b' integral and >= 2.
// 13 bin_packing: a <= form with beta' integral, >= 2 and equal to some |c_k|.
// 14 knapsack: a <= form with beta' integral and >= 2.
// 15 integer_knapsack: a <= form c^T x <= beta, all columns binary or integer with lb >= 0, all
//    c_k > 0 and integral, beta integral and >= 0 (no complementing).
// 16 mixed_binary: all columns binary or continuous (the table's form with no binary term, so
//    rows of continuous columns only count here too; they are also reported separately).
// 17 general_linear: the rest.
// Choices where the table leaves room: each row counts once (a ranged row may fit a class with
// either side), and >= sides are negated to their <= form for classes 8, 11 and 13-15, so
// x1 + x2 + x3 >= 2 is a set packing of the complements (at most one of them is 0).
enum class RowClass : uint8_t {
  Empty,
  Free,
  Singleton,
  Aggregation,
  Precedence,
  VariableBound,
  SetPartitioning,
  SetPacking,
  SetCovering,
  Cardinality,
  InvariantKnapsack,
  EquationKnapsack,
  BinPacking,
  Knapsack,
  IntegerKnapsack,
  MixedBinary,
  GeneralLinear,
};
inline constexpr int kNumRowClasses = 17;

const char* to_string(RowClass c);  // the JSON key: "empty", "set_partitioning", ...

std::vector<RowClass> classify_rows(const Model& model);

// The profile of a solve of `model` with `params` (as passed to solve(); the incumbent trace
// needs trace_incumbents) and its result; `version`: the solver version to record.
std::string profile_json(const Model& model, const SolverParams& params, const SolveResult& result,
                         const std::string& version);
void write_profile(const std::string& path, const Model& model, const SolverParams& params, const SolveResult& result,
                   const std::string& version);

}  // namespace milp
