#include "milp/presolve/clique_merge.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>

namespace milp {

namespace {

bool binary(const MipProblem& p, Index j) { return p.is_int[j] && p.col_lb[j] == 0.0 && p.col_ub[j] == 1.0; }

uint64_t edge_key(Index a, Index b) {
  if (a > b) std::swap(a, b);
  return (static_cast<uint64_t>(static_cast<uint32_t>(a)) << 32) | static_cast<uint32_t>(b);
}

}  // namespace

Index merge_cliques(MipProblem& p, Index max_new_rows) {
  const Index n = p.num_cols(), m = p.num_rows();
  const SparseMatrix R = p.A.transpose();
  // Set-packing rows: all binary, all coefficients 1, row_ub == 1, no lower side.
  std::vector<char> packing(m, 0);
  std::vector<std::vector<Index>> adj(n);
  std::unordered_map<uint64_t, Index> pair_row;  // edge -> pairwise row
  constexpr Index kMaxEdgeRow = 64;
  for (Index i = 0; i < m; ++i) {
    if (p.row_ub[i] != 1.0 || is_finite(p.row_lb[i]) && p.row_lb[i] > 0.0) continue;
    const Index len = R.end(i) - R.begin(i);
    if (len < 2 || len > kMaxEdgeRow) continue;
    bool ok = true;
    for (Index k = R.begin(i); k < R.end(i) && ok; ++k) ok = binary(p, R.index[k]) && R.value[k] == 1.0;
    if (!ok) continue;
    packing[i] = 1;
    for (Index a = R.begin(i); a < R.end(i); ++a)
      for (Index b = a + 1; b < R.end(i); ++b) {
        adj[R.index[a]].push_back(R.index[b]);
        adj[R.index[b]].push_back(R.index[a]);
      }
    if (len == 2) pair_row.emplace(edge_key(R.index[R.begin(i)], R.index[R.begin(i) + 1]), i);
  }
  if (pair_row.size() < 3) return 0;
  for (auto& a : adj) {
    std::sort(a.begin(), a.end());
    a.erase(std::unique(a.begin(), a.end()), a.end());
  }
  auto adjacent = [&](Index u, Index v) { return std::binary_search(adj[u].begin(), adj[u].end(), v); };

  // Grow a maximal clique from every uncovered pairwise row.
  std::vector<char> remove(m, 0);
  std::vector<std::vector<Index>> new_rows;
  for (const auto& [key, row] : pair_row) {
    if (remove[row] || static_cast<Index>(new_rows.size()) >= max_new_rows) continue;
    const Index u = static_cast<Index>(key >> 32), v = static_cast<Index>(key & 0xffffffffu);
    std::vector<Index> K{u, v};
    // Candidates: common neighbours of u and v, highest degree first.
    std::vector<Index> cand;
    std::set_intersection(adj[u].begin(), adj[u].end(), adj[v].begin(), adj[v].end(), std::back_inserter(cand));
    std::sort(cand.begin(), cand.end(), [&](Index a, Index b) { return adj[a].size() > adj[b].size(); });
    for (Index w : cand) {
      bool all = true;
      for (Index k : K) all = all && adjacent(k, w);
      if (all) K.push_back(w);
    }
    if (K.size() < 3) continue;
    // Pairwise rows inside K are dominated by the clique row.
    Index covered = 0;
    for (size_t a = 0; a < K.size(); ++a)
      for (size_t b = a + 1; b < K.size(); ++b) {
        auto it = pair_row.find(edge_key(K[a], K[b]));
        if (it != pair_row.end() && !remove[it->second]) remove[it->second] = 1, ++covered;
      }
    if (covered >= 2) {
      std::sort(K.begin(), K.end());
      new_rows.push_back(std::move(K));
    }
  }
  if (new_rows.empty()) return 0;

  // Rebuild the matrix row-wise: kept rows, then clique rows.
  std::vector<Index> J, I;
  std::vector<double> V, lb, ub;
  Index out = 0, removed = 0;
  for (Index i = 0; i < m; ++i) {
    if (remove[i]) {
      ++removed;
      continue;
    }
    for (Index k = R.begin(i); k < R.end(i); ++k) J.push_back(R.index[k]), I.push_back(out), V.push_back(R.value[k]);
    lb.push_back(p.row_lb[i]);
    ub.push_back(p.row_ub[i]);
    ++out;
  }
  for (const auto& K : new_rows) {
    for (Index j : K) J.push_back(j), I.push_back(out), V.push_back(1.0);
    lb.push_back(-kInf);
    ub.push_back(1.0);
    ++out;
  }
  p.A = SparseMatrix::from_triplets(n, out, J, I, V);
  p.row_lb = std::move(lb);
  p.row_ub = std::move(ub);
  return removed;
}

}  // namespace milp
