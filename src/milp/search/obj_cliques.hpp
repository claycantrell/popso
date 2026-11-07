#pragma once
// Objective propagation with clique implications (SCIP's pseudo-objective propagator with implications).
//
// Improving points satisfy c x <= U (the cutoff). Over the current bounds, P = sum_j min(c_j l_j, c_j u_j) is a
// lower bound on c x. Setting a binary literal to 1 forces every literal it conflicts with to 0, and each of those
// may move its column off its cheaper bound: if P plus the objective increase of the literal and of the literals
// it forces exceeds U, the literal is 0 in every improving point (its column is fixed to the other value). This is
// a reduction against the incumbent, like reduced-cost fixing: it keeps every point better than the cutoff, not
// the incumbent itself.
//
// Conflicts come from the rows' knapsack relaxations: a row side over binary literals with positive weights (a
// negative coefficient is the complemented literal), the other columns at their minimum activity bound; literals
// i and j conflict when w_i + w_j exceeds the slack. Sorted by weight, the literals that conflict with anything
// form a prefix (a list), and the partners of a literal in it are again a prefix of the list. So a row a x + C y
// <= C with C = sum a gives y = 1 => x = 0 for every item, which the prefix clique of the row (the separators'
// and the presolve's) only has for the heaviest item.
//
// hypothyroid-k1: SCIP's first root LP is -2856.9 instead of -2924.1 because of this propagation (with the
// incumbent -2692 of its vbounds heuristic before the root LP); on our presolved model, 67 fixings at that cutoff
// give LP -2851.08 and 74 at the optimum's cutoff give -2851.0 (root LP -2902.85 without).

#include <vector>

#include "milp/search/mip_problem.hpp"

namespace milp {

class ObjectiveCliques {
 public:
  // Conflict lists of the rows of p (at p's bounds), at most max_entries list entries in all (rows in order).
  ObjectiveCliques(const MipProblem& p, long long max_entries);

  Index lists() const { return static_cast<Index>(start_.size()) - 1; }
  long long entries() const { return static_cast<long long>(lit_.size()); }
  bool truncated() const { return truncated_; }

  // Fixes binaries in [lb, ub] (columns of p, lb <= ub) for the cutoff U on c x (no offset): changed columns are
  // appended to `changed`. Passes repeat (at most max_passes) while a pass fixed a column off its cheaper bound
  // (P grew). Only literals whose increase at p's bounds exceeds the room U - P are visited (in descending order
  // of it), so a call is cheap while the cutoff is far from the pseudo objective. work: columns and list entries
  // visited; no further literal once it reaches work_limit. Returns the number of fixings.
  Index propagate(double cutoff, std::vector<double>& lb, std::vector<double>& ub, std::vector<Index>& changed,
                  long long work_limit, int max_passes = 4, long long* work = nullptr);
  long long unions() const { return unions_; }  // literals that needed the exact union of their partners

 private:
  double literal_cost(Index l, bool own) const;
  Index partners(Index e) const;
  const std::vector<double>& cost_;
  const std::vector<char>& is_int_;
  Index n_ = 0;
  // List r: entries [start_[r], start_[r+1]) by weight descending; literal 2j (x_j = 1) or 2j+1 (x_j = 0).
  // Entries e, f of a list conflict iff w_[e] + w_[f] > slack_[r] (the tolerance is included in slack_).
  std::vector<Index> start_, lit_;
  std::vector<double> w_, slack_;
  std::vector<Index> list_of_;  // entry -> its list
  // Literal -> its entries (CSR over the 2n literals).
  std::vector<Index> mstart_, ment_;
  std::vector<double> ub0_;    // per literal: its increase at p's bounds, its lists added up
  std::vector<Index> order_;   // the literals in some list, by ub0_ descending
  bool truncated_ = false;
  long long unions_ = 0;
  // Scratch of propagate.
  std::vector<double> inc_;   // per literal: objective increase when forced to 0
  std::vector<char> true_;    // per literal: true at the pass's bounds
  std::vector<double> pre_;   // per list r: prefix sums of inc_ at start_[r] + r + k (k = 0..length)
  std::vector<int> ptrue_;    // ... and of true_
  std::vector<int> list_pass_;  // pass in which a list's prefix sums were built
  int pass_id_ = 0;
  std::vector<int> stamp_;
  std::vector<char> stamp_val_;
  int stamp_id_ = 0;
};

}  // namespace milp
