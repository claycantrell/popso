#pragma once
// Symmetry detection: generators of (a subgroup of) the group of column
// permutations that, together with a matching row permutation, map a MIP onto
// itself (same costs, bounds, types, row sides and coefficients).
//
// The problem is encoded as a colored bipartite graph (column vertices colored
// by cost/bounds/type, row vertices by their sides, edges by coefficient) and
// its automorphisms are searched with partition refinement plus
// individualization (nauty/saucy style, generators only, no canonical form).
// Every generator is verified exactly against the problem data before it is
// accepted, so hash collisions or budget limits can only lose symmetries.

#include <cstdint>
#include <utility>
#include <vector>

#include "milp/search/mip_problem.hpp"

namespace milp {

struct SymmetryOptions {
  long long work_limit = 0;   // elementary steps (0 = automatic, proportional to the problem size)
  double time_limit = 10.0;   // safety net in seconds (the work limit keeps runs deterministic)
  int max_generators = 1000;
  long long max_support = 5000000;  // total moved columns over all generators
};

struct Symmetry {
  // Generator g moves the columns gens[g][k].first to gens[g][k].second.
  std::vector<std::vector<std::pair<Index, Index>>> gens;
  // For each column: the generators moving it and its image (CSR by column).
  std::vector<Index> mov_start;
  std::vector<std::pair<int, Index>> movers;
  // Stabilizer chain: base[l] is the column individualized at level l of the
  // detection's first path; generator g fixes base[0..gen_level[g]-1], so the
  // generators of levels >= l generate (a subgroup of) the pointwise
  // stabilizer of base[0..l-1]. tree[l] is the Schreier tree of base[l] under
  // them (empty for levels without generators of their own); the product of
  // the tree sizes is the group order when the detection was complete.
  std::vector<Index> base;
  std::vector<int> gen_level;
  struct TreeEntry {
    Index point, parent;
    int gen;  // parent -> point (-1 at the root)
  };
  std::vector<std::vector<TreeEntry>> tree;
  double seconds = 0.0;
  long long work = 0;
  bool aborted = false;   // budget exhausted (the generators found are still valid)
  Index moved_cols = 0;   // columns moved by some generator

  bool empty() const { return gens.empty(); }
  int num_generators() const { return static_cast<int>(gens.size()); }
  double log10_order() const;
  // Dense t in the group of the generators of levels >= l with t(base[l]) = y;
  // false if y is not in that orbit.
  bool transversal(int l, Index y, std::vector<Index>& t) const;
  void build_chain(Index n);
};

// Scratch space for PermGroup::stabilizer (one per thread).
struct GroupScratch {
  std::vector<std::vector<Index>> dense;
  std::vector<Index> parent, ux, uz, uzinv, g;
  std::vector<int> via;
  std::vector<Index> orbit;
};

// A group of column permutations given by generators, each a list of
// (column, image) pairs for the moved columns, sorted by column.
struct PermGroup {
  std::vector<std::vector<std::pair<Index, Index>>> gens;

  Index image(int g, Index j) const;
  // Orbit of j (j first). `mark` must be all zero with one entry per column;
  // it is restored before returning.
  void orbit(Index j, std::vector<Index>& out, std::vector<char>& mark) const;
  // Generators of a subgroup of the stabilizer of j: the generators fixing j
  // plus random Schreier generators u_{s(x)}^-1 s u_x of the orbit of j (all
  // exact products of group elements, so every one fixes j and is in the
  // group). Deterministic for a given seed. n = number of columns.
  PermGroup stabilizer(Index j, Index n, uint64_t seed, GroupScratch& s, int max_kept = 12, int max_new = 12) const;
  long long bytes() const;
};

Symmetry detect_symmetry(const MipProblem& p, const SymmetryOptions& opts = {});

}  // namespace milp
