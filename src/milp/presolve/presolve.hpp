#pragma once
// MIP presolve: transforms a MipProblem into a smaller equivalent one and
// maps solutions of the reduced problem back (primal postsolve only).
//
// Reductions (all valid for MILP, applied in rounds until little progress):
//   - empty rows / empty columns, fixed columns
//   - singleton rows -> column bounds
//   - activity-based bound tightening, redundant and forcing rows
//   - coefficient strengthening on one-sided rows (integer columns)
//   - parallel rows; parallel columns with proportional costs (merged)
//   - dual fixing (lock-based); dominated columns
//   - doubleton equations (aggregation; integer pairs when integrality is preserved)
//   - substitution of (implied) free columns and column singletons out of equations
//   - probing on binaries: fixings, bound tightening, equivalences x = y / x = 1 - y, lifted rows
//     (repeated passes)
//   - sparsification: equations added to other rows to cancel entries
// Afterwards implied integers are detected and flagged in reduced().implied_int (their bounds rounded).

#include <atomic>
#include <utility>
#include <vector>

#include "milp/core/log.hpp"
#include "milp/core/types.hpp"
#include "milp/search/mip_problem.hpp"

namespace milp {

struct PresolveOptions {
  int max_rounds = 25;
  double feas_tol = 1e-6;
  bool bound_tightening = true;
  bool coef_tightening = true;
  bool parallel_rows = true;
  bool parallel_cols = true;
  bool dual_fixing = true;
  bool dominated_cols = true;
  bool aggregation = true;
  bool substitution = true;  // free column / column singleton substitution
  // Also substitute non-singleton continuous columns with a nonzero objective coefficient whose bounds the
  // equation does not imply (the bounds stay as the row's sides; MILP_PRESOLVE_ON=subst_bounded). Off: such a
  // column keeps its bounds and variable bounds for the cut separators (tr12-30: the stock columns of the
  // inventory balances, with holding costs, became row slacks and the c-MIR path aggregation lost the
  // (l,S)-like cuts; root 127730 -> 130056 without the substitution). Zero-cost columns are still substituted
  // (empirical: neos-3627168-kasai's stocks, cost 0, solve 2-3x faster substituted).
  bool substitution_bounded = false;
  bool probing = true;
  // Lifting of one-sided rows by probing implications (probing.hpp; MILP_PRESOLVE_OFF=probing_lift).
  bool probing_lift = true;
  // Probing passes: another one (half the work limit of the previous) while a pass fixes, substitutes or
  // lifts something (MILP_PRESOLVE_OFF=probing_passes: a single pass; MILP_PROBING_PASSES=<n>).
  int probing_passes = 3;
  // Implications between binaries found by probing are kept and propagated by later probes and passes
  // (MILP_PRESOLVE_OFF=probing_implications). Opt-in (MILP_PRESOLVE_ON=probing_cliques): those of the form
  // x = 1 => y = 0 join the conflict graph of the clique extension (puriri: LP 60353 without, 60337 with, and
  // 72k vs 81k nonzeros; mzzv11 35.8% vs 35.5% of the root gap closed, 112k vs 130k nonzeros: the liftings
  // already add the members that matter).
  bool probing_implications = true;
  bool probing_cliques = false;
  // No substitution of binaries out of equations before the first probing pass (MILP_PRESOLVE_OFF=probing_keep).
  bool probing_keep_binaries = true;
  // Sparsification: equations added to rows that share columns with them when that cancels more entries than it
  // adds (PresolveWork::sparsify; MILP_PRESOLVE_OFF=sparsify).
  bool sparsify = true;
  // ... without equations that contain a continuous column with a variable bound (an inequality with continuous
  // columns of one sign and integer columns of the other): the cut separators aggregate through them
  // (MILP_PRESOLVE_OFF=sparsify_protect).
  bool sparsify_protect = true;
  bool clique_merging = true;
  bool clique_extension = true;  // extend set-packing rows by the conflict graph of all rows
  bool implied_int = true;       // implied-integer detection on the reduced problem (presolve/implied_int.hpp)
  bool implied_int_dual = true;  // ... including the dual rule
  double probing_effort = 1.0;  // scales the probing work limit
  // Optional (sub-MIPs): no further round or step once *stop is set or time_limit seconds have passed since the
  // start; the reductions so far are kept. The top-level presolve sets neither (work limits only).
  const std::atomic<bool>* stop = nullptr;
  double time_limit = kInf;

  // Defaults with switches from MILP_PRESOLVE_OFF, a comma-separated list of option names to disable
  // (e.g. "probing,dominated"); "all" disables every reduction except the basic ones.
  static PresolveOptions from_env();
};

enum class PresolveStatus {
  Reduced,     // reduced() holds an equivalent (possibly empty) problem
  Infeasible,  // the problem was proven infeasible
  Unbounded,   // unbounded if feasible (i.e. infeasible or unbounded)
};

struct PresolveStats {
  Index rows_before = 0, cols_before = 0, nnz_before = 0;
  Index rows_after = 0, cols_after = 0, nnz_after = 0;
  long long bound_changes = 0;
  long long coef_changes = 0;
  long long fixed_cols = 0;
  long long removed_rows = 0;
  long long parallel_rows = 0;
  long long parallel_cols = 0;
  long long dominated = 0;
  long long dual_bounds = 0;  // dual bound strengthening
  long long aggregations = 0;
  long long substitutions = 0;
  long long probes = 0, probing_fixed = 0, probing_bounds = 0, probing_equivs = 0;
  long long probing_lifts = 0, probing_lift_fill = 0;  // lifted rows; of those, new entries
  long long probing_conflicts = 0;                     // implications between binaries kept (literal pairs)
  int probing_passes = 0;
  long long sparsify_rows = 0, sparsify_removed = 0;  // row combinations; nonzeros removed by them (net)
  long long sparsify_protected = 0;                    // equations not used for their variable-bound columns
  int sparsify_runs = 0;
  long long clique_extended = 0, clique_dominated = 0, clique_fixed = 0;
  long long implied_primal = 0, implied_dual = 0, implied_int = 0;  // new by rule; flagged in total
  int rounds = 0;
  double time = 0.0;
};

class Presolver {
 public:
  PresolveStatus run(const MipProblem& in, const PresolveOptions& opt, const Logger& log);

  const MipProblem& reduced() const { return reduced_; }
  // Maps a solution of reduced() to the original column space.
  std::vector<double> postsolve(const std::vector<double>& x_reduced) const;
  // Maps a point of the original column space to reduced() (inverse of postsolve on the kept
  // columns). The result need not be feasible: dual reductions may have cut the point off.
  std::vector<double> crush(const std::vector<double>& x) const;

  const PresolveStats& stats() const { return stats_; }
  // Original column index of each reduced column.
  const std::vector<Index>& col_map() const { return col_map_; }

  // One postsolve step. Fixed: x[col] = rhs. Substitution: x[col] = (rhs - sum terms) / coef, clamped to
  // [lb, ub] (the column's bounds when it was substituted; a clamp only moves rounding errors into the row).
  //
  // Split (merged parallel columns, other >= 0): x[other] holds y = x[other] + coef * x[col]; it is split into
  // x[col] in [lb, ub] and x[other] in [olb, oub] (integral values if integer).
  struct Step {
    Index col = -1;
    double rhs = 0.0, coef = 1.0;
    std::vector<std::pair<Index, double>> terms;  // empty for fixed columns
    double lb = -kInf, ub = kInf;
    Index other = -1;
    double olb = -kInf, oub = kInf;
    bool integer = false;
  };

 private:
  friend class PresolveWork;
  MipProblem reduced_;
  PresolveStats stats_;
  Index orig_cols_ = 0;
  std::vector<Index> col_map_;
  std::vector<char> orig_is_int_;
  std::vector<Step> stack_;  // applied in reverse order by postsolve
};

}  // namespace milp
