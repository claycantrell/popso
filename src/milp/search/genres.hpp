#pragma once
// Cut-based (generalized-resolution) conflict analysis, after Mexi, Serrano,
// Berthold, Gleixner, Nordstrom, "Cut-based Conflict Analysis in Mixed Integer
// Programming", INFORMS J. Computing (2025), doi:10.1287/ijoc.2024.0999.
//
// Instead of a clause over bound literals (analyze_propagation_conflict), the
// learned object is a linear constraint: starting from the constraint that
// failed (a row side, a conflict, or the LP's dual proof), the latest bound
// change it needs is eliminated by a nonnegative combination with that bound's
// reason, repeatedly, until the constraint is asserting (one needed bound at
// the highest decision level: first UIP) or the needed bound is a decision.
// Each step keeps the constraint infeasible under the bounds still in effect;
// when the plain combination does not (the reason propagated an integer bound
// after rounding), the reason is reduced first: continuous bounds it relies on
// are resolved out with their own reasons, the rest are weakened to their
// global bounds, and a c-MIR rounding with the resolved variable's coefficient
// as divisor (variables complemented at the global bound their local bound
// sits at) makes it propagate tightly. Every step is a valid derivation from
// the rows, conflicts and global bounds, so the result holds for every
// solution better than the cutoff whatever the trail says; the trail only
// guides the derivation.

#include <vector>

#include "milp/core/sparse.hpp"
#include "milp/core/types.hpp"
#include "milp/search/conflict.hpp"
#include "milp/search/trail.hpp"

namespace milp {

struct GenresOptions {
  Index max_nnz = 1 << 30;  // the learned constraint (and intermediate ones up to 2x) may not be denser
  int max_steps = 100;      // resolution steps
  long long max_work = 200000;  // nonzeros touched
  bool partial = true;      // on failure after >= 1 step, return the last (valid, infeasible) resolvent
  bool coef_tightening = true;
};

struct GenresInfo {
  enum Outcome { kUip, kDecision, kGlobal, kPartial, kFail };
  enum Reason {
    kNone,
    kNoStart,      // no usable failing constraint (clause on non-binaries, infinite bounds, not infeasible)
    kNoStep,       // the first needed bound is a decision: nothing to resolve
    kClauseReason, // a reason is a clause on general integers
    kReduction,    // the (reduced) resolvent is not infeasible
    kContDecision, // a continuous bound the reason relies on is a decision
    kLimit,        // nonzeros, steps or work exceeded
    kNumerics,
    kFinalCheck    // the cleaned result does not prove the node infeasible
  };
  Outcome outcome = kFail;
  Reason reason = kNone;
  int steps = 0;      // resolution steps
  int reduced = 0;    // steps that needed a reduction of the reason
  int mir = 0;        // of these, steps with a c-MIR rounding
};

// Analyzes the failure recorded in `trail` (rows: the propagator's rows with
// sides [row_lb, row_ub]; reasons of kConflict entries are the conflicts they
// point to). Writes a linear conflict to `out` (cutoff: the smallest of the
// conflicts used; the caller must also apply the cutoff that rows such as the
// objective row and the global bounds [glb, gub] rely on). Returns true if a
// constraint was learned (outcome kUip, kDecision, kGlobal or kPartial).
bool analyze_genres(const BoundTrail& trail, const SparseMatrix& rows, const std::vector<double>& row_lb,
                    const std::vector<double>& row_ub, const std::vector<char>& is_int, const std::vector<double>& glb,
                    const std::vector<double>& gub, const GenresOptions& opt, Conflict& out, GenresInfo* info = nullptr);

}  // namespace milp
