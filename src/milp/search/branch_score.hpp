#pragma once
// Branching scores and the strong-branching stopping rule.
//
// Score corrections after P. Shah, S. S. Dey, "Improving Full Strong Branching Decisions by
// Incorporating Additional Information", arXiv:2507.09455 (2025):
//  - efficacious gains: a child's objective gain counts at most the primal-dual gap G of the node
//    (cutoff - node LP objective): a gain beyond G prunes the child, more does not help; an
//    infeasible child counts G;
//  - min/max exponents: q_min^0.3 q_max^0.7 instead of the product q_min q_max;
//  - last-assignment asymmetry: with a = (n_dn - n_up) / (n_dn + n_up) over the leaves whose last
//    branching went down / up, the score gets a factor q_up^(eta a) if a > 0, else q_dn^(-eta a)
//    (eta = 0.15); the leaves counted are the infeasible or integral ones, or all leaves while
//    those are under 5% of all (the paper's PA-LA rule).
//
// Probabilistic lookahead after G. Mexi, S. Shamsi, M. Besancon, P. Le Bodic, "Probabilistic
// Lookahead Strong Branching via a Stochastic Abstract Branching Model", CPAIOR 2024
// (arXiv:2312.07041): the geometric-mean gains g of the candidates probed at a node are a sample of
// a distribution with a point mass p0 at zero and a fitted continuous part (exponential or Pareto,
// maximum likelihood); a candidate with gain g closes the gap G with a single-variable tree of depth
// d = ceil(G / g), 2^(d+1) - 1 nodes, and each strong-branching call costs 2 nodes. Probing one more
// candidate is worth it iff the expected tree size afterwards (best depth min(d*, D) of the new
// sample D, plus its 2 nodes) is below the current one: sum_{d < d*} P(D = d) (2^d* - 2^d) > 1.

#include <algorithm>
#include <cmath>

#include "milp/core/types.hpp"

namespace milp {

struct BranchScoreRule {
  bool cap = false;     // gains capped at gap (efficacious gains)
  bool minmax = false;  // q_min^0.3 q_max^0.7 instead of the product
  double gap = kInf;    // cutoff - node objective (kInf: no incumbent, nothing is capped)
  double a_dn = 0.0, a_up = 0.0;  // asymmetry exponents (at most one positive)
  static constexpr double kEps = 1e-6;
  bool plain() const { return !cap && !minmax && a_dn == 0.0 && a_up == 0.0; }
};

// Score of a candidate with (estimated or probed) child gains dn, up >= 0 (kInf: child pruned).
// The plain rule is the product max(dn, eps) max(up, eps).
inline double branch_score(double dn, double up, const BranchScoreRule& r) {
  constexpr double eps = BranchScoreRule::kEps;
  double qd = std::max(dn, eps), qu = std::max(up, eps);
  if (r.plain()) return qd * qu;
  if (r.cap) {
    const double cap = std::max(r.gap, eps);
    qd = std::min(qd, cap), qu = std::min(qu, cap);
  }
  double s = r.minmax ? std::pow(std::min(qd, qu), 0.3) * std::pow(std::max(qd, qu), 0.7) : qd * qu;
  if (r.a_dn > 0.0) s *= std::pow(qd, r.a_dn);
  if (r.a_up > 0.0) s *= std::pow(qu, r.a_up);
  return s;
}

// Last-assignment exponents from the leaf counts by the direction of their last branching:
// inf_* = infeasible or integral leaves, all_* = all leaves (pruned by bound too). PA-LA: the
// infeasible / integral leaves decide once there are at least min_leaves of them and they make up
// at least 5% of all leaves; else all leaves (rebalancing), if there are min_leaves of those.
inline void asymmetry_exponents(double inf_dn, double inf_up, double all_dn, double all_up, double eta,
                                double min_leaves, double* a_dn, double* a_up) {
  *a_dn = *a_up = 0.0;
  const double ni = inf_dn + inf_up, na = all_dn + all_up;
  double dn, up;
  if (ni >= min_leaves && ni >= 0.05 * na) dn = inf_dn, up = inf_up;
  else if (na >= min_leaves) dn = all_dn, up = all_up;
  else return;
  const double a = (dn - up) / (dn + up);
  if (a > 0.0) *a_up = eta * a;
  else *a_dn = -eta * a;
}

// Probabilistic lookahead: the gains of the candidates probed at one node.
class LookaheadSample {
 public:
  void clear() { k_ = zeros_ = 0, sum_ = sum_log_ = 0.0, min_ = kInf; }
  // Geometric-mean gain of a probed candidate (both children finite).
  void add(double g) {
    if (!(g > 1e-9)) {
      ++zeros_;
      return;
    }
    ++k_;
    sum_ += g;
    sum_log_ += std::log(g);
    min_ = std::min(min_, g);
  }
  int nonzero() const { return k_; }
  int size() const { return k_ + zeros_; }

  // P(g >= x) under the fitted distribution (dist 0: exponential, 1: Pareto).
  double survival(double x, int dist) const {
    const int n = k_ + zeros_;
    if (n == 0 || k_ == 0) return 0.0;
    const double nz = static_cast<double>(k_) / n;
    if (x <= 0.0) return nz;
    if (dist == 0) return nz * std::exp(-x * k_ / sum_);
    if (x <= min_) return nz;
    // MLE: x_m = min, alpha = k / sum ln(g / x_m); all gains equal: a point mass at x_m.
    const double s = sum_log_ - k_ * std::log(min_);
    if (s <= 1e-12) return 0.0;
    const double alpha = k_ / s;
    return nz * std::exp(-alpha * std::log(x / min_));
  }

  // Whether probing one more candidate is expected to shrink the tree (see the file comment), for
  // the gap G and the best geometric-mean gain best_g found so far.
  bool worth_more(double G, double best_g, int dist) const {
    if (!(G > 0.0) || !std::isfinite(G)) return true;  // no incumbent: no model
    if (best_g >= G) return false;                      // the best candidate closes the gap
    if (!(best_g > 0.0)) return survival(0.0, dist) > 0.0 || size() == 0;
    const double ds = std::ceil(G / best_g);
    if (ds > 1000.0) return survival(G / (ds - 1.0), dist) > 0.0;
    const int dstar = static_cast<int>(ds);
    if (dstar <= 1) return false;
    // V = sum_{d < d*} P(D = d) (1 - 2^(d - d*)), P(D <= d) = P(g >= G / d).
    const double p_less = survival(G / (dstar - 1), dist);
    double near = 0.0;
    for (int d = std::max(1, dstar - 60); d < dstar; ++d) {
      const double pd = survival(G / d, dist) - (d > 1 ? survival(G / (d - 1), dist) : 0.0);
      near += std::max(pd, 0.0) * std::ldexp(1.0, d - dstar);
    }
    return p_less - near > std::ldexp(1.0, -dstar);
  }

 private:
  int k_ = 0, zeros_ = 0;
  double sum_ = 0.0, sum_log_ = 0.0, min_ = kInf;
};

}  // namespace milp
