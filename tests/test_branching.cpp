// Tests for the branching score corrections (Shah-Dey) and the probabilistic lookahead stopping
// rule (Mexi et al.) in search/branch_score.hpp.

#include <cmath>
#include <vector>

#include "milp/search/branch_score.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

// Expected tree size of the abstract model after one more strong-branching call minus the current
// one, computed directly: t_i = 2^(d*+1) - 1 + 2i, E[t_{i+1}] = sum_d (2^(d+1) - 1) P(min(d*, D) = d)
// + 2(i+1) with P(D <= d) = P(g >= G / d).
double expected_change(const LookaheadSample& s, double G, double best_g, int dist) {
  const int dstar = static_cast<int>(std::ceil(G / best_g));
  double e = 0.0, below = 0.0;
  for (int d = 1; d < dstar; ++d) {
    const double pd = s.survival(G / d, dist) - (d > 1 ? s.survival(G / (d - 1), dist) : 0.0);
    e += (std::ldexp(1.0, d + 1) - 1.0) * pd;
    below += pd;
  }
  e += (std::ldexp(1.0, dstar + 1) - 1.0) * (1.0 - below);
  return e + 2.0 - (std::ldexp(1.0, dstar + 1) - 1.0);
}

}  // namespace

TEST(branch_score_plain_is_product) {
  BranchScoreRule r;
  CHECK_NEAR(branch_score(2.0, 3.0, r), 6.0, 1e-12);
  CHECK_NEAR(branch_score(0.0, 3.0, r), 3e-6, 1e-15);
  // Plain scores do not cap: a gain beyond any gap counts in full.
  r.gap = 1.0;
  CHECK_NEAR(branch_score(10.0, 3.0, r), 30.0, 1e-12);
}

TEST(branch_score_capped_gains) {
  BranchScoreRule r;
  r.cap = true;
  r.gap = 5.0;
  CHECK_NEAR(branch_score(10.0, 3.0, r), 15.0, 1e-12);
  CHECK_NEAR(branch_score(kInf, 3.0, r), 15.0, 1e-12);  // a pruned child counts the gap
  CHECK_NEAR(branch_score(kInf, kInf, r), 25.0, 1e-12);
  // Without an incumbent nothing is capped.
  r.gap = kInf;
  CHECK_NEAR(branch_score(10.0, 3.0, r), 30.0, 1e-12);
  // A capped candidate no longer beats a balanced one it would win against in the product.
  r.gap = 4.0;
  BranchScoreRule plain;
  CHECK(branch_score(100.0, 1.0, plain) > branch_score(3.0, 3.0, plain));
  CHECK(branch_score(100.0, 1.0, r) < branch_score(3.0, 3.0, r));
}

TEST(branch_score_minmax_exponents) {
  BranchScoreRule r;
  r.minmax = true;
  CHECK_NEAR(branch_score(2.0, 8.0, r), std::pow(2.0, 0.3) * std::pow(8.0, 0.7), 1e-12);
  CHECK_NEAR(branch_score(8.0, 2.0, r), std::pow(2.0, 0.3) * std::pow(8.0, 0.7), 1e-12);
  r.cap = true;
  r.gap = 4.0;
  CHECK_NEAR(branch_score(2.0, 8.0, r), std::pow(2.0, 0.3) * std::pow(4.0, 0.7), 1e-12);
}

TEST(branch_score_asymmetry) {
  double a_dn = -1, a_up = -1;
  // Too few leaves: no correction.
  asymmetry_exponents(2, 3, 4, 4, 0.15, 10, &a_dn, &a_up);
  CHECK(a_dn == 0.0 && a_up == 0.0);
  // Infeasible / integral leaves decide (>= 10 and >= 5% of all): more after down branchings.
  asymmetry_exponents(30, 10, 200, 200, 0.15, 10, &a_dn, &a_up);
  CHECK_NEAR(a_up, 0.15 * 0.5, 1e-12);
  CHECK(a_dn == 0.0);
  // Under 5% of all leaves: all leaves decide (rebalancing), here more after up branchings.
  asymmetry_exponents(6, 6, 100, 300, 0.15, 10, &a_dn, &a_up);
  CHECK_NEAR(a_dn, 0.15 * 0.5, 1e-12);
  CHECK(a_up == 0.0);
  BranchScoreRule r;
  r.a_up = 0.1;
  CHECK_NEAR(branch_score(2.0, 3.0, r), 6.0 * std::pow(3.0, 0.1), 1e-12);
}

TEST(lookahead_fit_and_survival) {
  LookaheadSample s;
  CHECK(s.survival(1.0, 0) == 0.0);
  for (double g : {1.0, 2.0, 3.0}) s.add(g);
  s.add(0.0);
  CHECK(s.nonzero() == 3 && s.size() == 4);
  // Exponential: point mass 1/4 at zero, rate 1 / mean = 0.5.
  CHECK_NEAR(s.survival(0.0, 0), 0.75, 1e-12);
  CHECK_NEAR(s.survival(2.0, 0), 0.75 * std::exp(-1.0), 1e-12);
  // Pareto: x_m = 1, alpha = 3 / (ln 2 + ln 3).
  const double alpha = 3.0 / std::log(6.0);
  CHECK_NEAR(s.survival(0.5, 1), 0.75, 1e-12);
  CHECK_NEAR(s.survival(4.0, 1), 0.75 * std::pow(4.0, -alpha), 1e-12);
}

TEST(lookahead_stopping_rule) {
  LookaheadSample s;
  for (double g : {1.0, 1.5, 0.5, 2.0, 1.0}) s.add(g);
  // No incumbent (no gap): no model, keep probing; a candidate that closes the gap: stop.
  CHECK(s.worth_more(kInf, 1.0, 0));
  CHECK(!s.worth_more(1.5, 2.0, 0));
  // Far from closing the gap: another candidate may shorten the tree a lot.
  CHECK(s.worth_more(20.0, 2.0, 0));
  // Gains far below what would change d*: stop (d* = 2 needs g >= 10 to improve).
  LookaheadSample small;
  for (double g : {0.01, 0.02, 0.015, 0.01}) small.add(g);
  CHECK(!small.worth_more(10.0, 5.0, 0));
  // The rule agrees with the direct expected tree-size comparison.
  for (int dist = 0; dist < 2; ++dist) {
    for (double G : {3.0, 7.5, 12.0, 40.0}) {
      for (double best : {0.8, 1.7, 2.5, 6.0}) {
        if (best >= G) continue;
        const bool direct = expected_change(s, G, best, dist) < 0.0;
        CHECK(s.worth_more(G, best, dist) == direct);
      }
    }
  }
}
