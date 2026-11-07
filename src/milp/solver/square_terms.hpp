#pragma once
// Separable convex quadratic objectives on integer columns as a MILP (solver.cpp, solve_quadratic).
//
// When Q is diagonal and every column with Q_jj != 0 is integer, each term 1/2 Q_jj x_j^2 has an exact
// linear representation: a binary has x^2 = x; for a general integer, with a center c and d = x - c,
//   d^2 = max over integers k of (2k + 1) d - k(k + 1)   at every integer d
// (the tangent at k is exact at d = k and d = k + 1). A column t_j >= the tangents for the offsets k in a
// set K_j stands for d_j^2, so 1/2 Q_jj x_j^2 = 1/2 Q_jj (t_j + 2 c x_j - c^2). With finitely many
// tangents the MILP is a relaxation of the model (each tangent lies below d^2), and it is exact at a
// solution whose d_j has the tangent k = d_j or k = d_j - 1 in K_j; solve_quadratic adds windows of
// tangents where a solution lacks them and solves again. The tangents give the convex hull of the
// parabola's integer points, a much tighter relaxation than the continuous x^2 of the MIQP search.

#include <set>
#include <vector>

#include "milp/core/model.hpp"

namespace milp {

class SquareTerms {
 public:
  // False unless Q is diagonal, every column with a nonzero entry is integer (or binary, or continuous but
  // implied integer: +-1 in an equality row with an integer right-hand side over integer columns with integer
  // coefficients; integer in the MILP) and the objective is convex in the model's sense; the model must have
  // no logical constraints. window: the tangent offsets
  // -window..window around a center (and a coarse grid of +-2^i beyond them).
  bool init(const Model& m, int window);
  // Centers at the rounded values of x (a relaxation's point; clipped to the bounds), with their windows.
  void center_at(const std::vector<double>& x);
  // Adds a window around d_j for every term whose tangents miss x's value; returns how many did.
  int refine(const std::vector<double>& x);
  bool exact_at(const std::vector<double>& x) const;
  // The MILP: the model's columns first (same order), then one t column per general-integer term.
  Model build() const;
  // A model point extended to the MILP's columns (t_j = d_j^2, its exact value).
  std::vector<double> extend(const std::vector<double>& x) const;
  size_t num_terms() const { return terms_.size(); }
  size_t num_general() const;
  size_t num_tangents() const;

 private:
  struct Term {
    Index col = -1;
    double half_q = 0.0;  // 1/2 Q_jj, model sense
    bool binary = false;
    long long center = 0;
    long long lo = 0, hi = 0;  // the offsets' range: d in [lo, hi] (from the column bounds, capped)
    std::set<long long> ks;    // tangent offsets
  };
  void add_window(Term& t, long long d);
  const Model* m_ = nullptr;
  std::vector<Term> terms_;
  std::vector<Index> implied_;  // continuous columns of terms, implied integer
  int window_ = 25;
};

}  // namespace milp
