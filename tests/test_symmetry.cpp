// Tests for symmetry detection (every generator must map the problem exactly
// onto itself; group orders of small graphs) and for orbital branching
// end to end (random symmetric models solve to the same optimum with and
// without symmetry handling).

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <random>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "milp/core/model.hpp"
#include "milp/search/mip_problem.hpp"
#include "milp/search/orbitope.hpp"
#include "milp/search/sst.hpp"
#include "milp/search/symmetry.hpp"
#include "milp/solver/solver.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

int scaled(int base) {
#ifdef _MSC_VER
#pragma warning(suppress : 4996)
#endif
  const char* v = std::getenv("MILP_TEST_SCALE");
  return v ? base * std::max(1, std::atoi(v)) : base;
}

void set_env(const char* name, bool on) {
#ifdef _WIN32
  _putenv_s(name, on ? "1" : "");
#else
  if (on) setenv(name, "1", 1);
  else unsetenv(name);
#endif
}

void set_env_value(const char* name, const char* value) {
#ifdef _WIN32
  _putenv_s(name, value);
#else
  setenv(name, value, 1);
#endif
}

using RowKey = std::tuple<double, double, std::vector<std::pair<Index, double>>>;

// Independent check: sigma maps columns onto columns with equal data, and the
// multiset of rows (sides + entries) is invariant under sigma.
bool maps_onto_itself(const MipProblem& p, const std::vector<Index>& sigma) {
  const Index n = p.num_cols(), m = p.num_rows();
  std::vector<char> seen(n, 0);
  for (Index j = 0; j < n; ++j) {
    const Index k = sigma[j];
    if (k < 0 || k >= n || seen[k]) return false;
    seen[k] = 1;
    if (p.cost[j] != p.cost[k] || p.col_lb[j] != p.col_lb[k] || p.col_ub[j] != p.col_ub[k] ||
        p.is_int[j] != p.is_int[k]) {
      return false;
    }
  }
  std::vector<RowKey> orig(m), img(m);
  for (Index r = 0; r < m; ++r) {
    std::get<0>(orig[r]) = std::get<0>(img[r]) = p.row_lb[r];
    std::get<1>(orig[r]) = std::get<1>(img[r]) = p.row_ub[r];
  }
  for (Index j = 0; j < n; ++j) {
    for (Index k = p.A.begin(j); k < p.A.end(j); ++k) {
      std::get<2>(orig[p.A.index[k]]).emplace_back(j, p.A.value[k]);
      std::get<2>(img[p.A.index[k]]).emplace_back(sigma[j], p.A.value[k]);
    }
  }
  for (Index r = 0; r < m; ++r) {
    std::sort(std::get<2>(orig[r]).begin(), std::get<2>(orig[r]).end());
    std::sort(std::get<2>(img[r]).begin(), std::get<2>(img[r]).end());
  }
  std::sort(orig.begin(), orig.end());
  std::sort(img.begin(), img.end());
  return orig == img;
}

std::vector<Index> generator(const Symmetry& s, int g, Index n) {
  std::vector<Index> sigma(n);
  for (Index j = 0; j < n; ++j) sigma[j] = j;
  for (const auto& [j, k] : s.gens[g]) sigma[j] = k;
  return sigma;
}

bool all_generators_valid(const MipProblem& p, const Symmetry& s) {
  for (int g = 0; g < s.num_generators(); ++g) {
    if (s.gens[g].empty() || !maps_onto_itself(p, generator(s, g, p.num_cols()))) return false;
  }
  return true;
}

// Group order from the stabilizer chain (exact if the detection was complete).
double group_order(const Symmetry& s) {
  double o = 1.0;
  for (const auto& t : s.tree) o *= t.empty() ? 1.0 : static_cast<double>(t.size());
  return o;
}

std::vector<Index> full_orbit(const Symmetry& s, Index j, Index n) {
  PermGroup g;
  g.gens = s.gens;
  std::vector<Index> out;
  std::vector<char> mark(n, 0);
  g.orbit(j, out, mark);
  std::sort(out.begin(), out.end());
  return out;
}

// Maximum independent set (as minimization) on a graph, optional weights.
Model independent_set(int nv, const std::vector<std::pair<int, int>>& edges, const std::vector<double>& w = {}) {
  ModelBuilder b;
  for (int v = 0; v < nv; ++v) b.add_col("x" + std::to_string(v), w.empty() ? -1.0 : -w[v], 0, 1, VarType::Binary);
  for (size_t e = 0; e < edges.size(); ++e) {
    const Index r = b.add_row("e" + std::to_string(e), -kInf, 1);
    b.add_coef(r, edges[e].first, 1);
    b.add_coef(r, edges[e].second, 1);
  }
  return b.build();
}

std::vector<std::pair<int, int>> cycle(int n, int off = 0) {
  std::vector<std::pair<int, int>> e;
  for (int i = 0; i < n; ++i) e.emplace_back(off + i, off + (i + 1) % n);
  return e;
}

// Bin packing with k identical bins.
Model bin_packing(const std::vector<int>& sizes, int k, int cap) {
  ModelBuilder b;
  const int n = static_cast<int>(sizes.size());
  std::vector<Index> y(k);
  std::vector<std::vector<Index>> x(n, std::vector<Index>(k));
  for (int c = 0; c < k; ++c) y[c] = b.add_col("y" + std::to_string(c), 1, 0, 1, VarType::Binary);
  for (int i = 0; i < n; ++i) {
    for (int c = 0; c < k; ++c) {
      x[i][c] = b.add_col("x" + std::to_string(i) + "_" + std::to_string(c), 0, 0, 1, VarType::Binary);
    }
  }
  for (int i = 0; i < n; ++i) {
    const Index r = b.add_row("a" + std::to_string(i), 1, 1);
    for (int c = 0; c < k; ++c) b.add_coef(r, x[i][c], 1);
  }
  for (int c = 0; c < k; ++c) {
    const Index r = b.add_row("c" + std::to_string(c), -kInf, 0);
    for (int i = 0; i < n; ++i) b.add_coef(r, x[i][c], sizes[i]);
    b.add_coef(r, y[c], -cap);
  }
  return b.build();
}

// Graph coloring with k identical colors: x[v][c], w[c]; min sum w.
Model coloring(int nv, const std::vector<std::pair<int, int>>& edges, int k) {
  ModelBuilder b;
  std::vector<Index> w(k);
  std::vector<std::vector<Index>> x(nv, std::vector<Index>(k));
  for (int c = 0; c < k; ++c) w[c] = b.add_col("w" + std::to_string(c), 1, 0, 1, VarType::Binary);
  for (int v = 0; v < nv; ++v) {
    for (int c = 0; c < k; ++c) {
      x[v][c] = b.add_col("x" + std::to_string(v) + "_" + std::to_string(c), 0, 0, 1, VarType::Binary);
    }
  }
  for (int v = 0; v < nv; ++v) {
    const Index r = b.add_row("a" + std::to_string(v), 1, 1);
    for (int c = 0; c < k; ++c) b.add_coef(r, x[v][c], 1);
  }
  for (size_t e = 0; e < edges.size(); ++e) {
    for (int c = 0; c < k; ++c) {
      const Index r = b.add_row("e" + std::to_string(e) + "_" + std::to_string(c), -kInf, 0);
      b.add_coef(r, x[edges[e].first][c], 1);
      b.add_coef(r, x[edges[e].second][c], 1);
      b.add_coef(r, w[c], -1);
    }
  }
  return b.build();
}

// k copies of a random block (integers in [0, ub], a couple of continuous
// columns) plus linking rows with equal coefficients on corresponding
// columns of every copy; columns and rows are emitted in a shuffled order.
Model duplicated_blocks(std::mt19937& rng, int k, int ni, int nc, int nr, int nlink) {
  std::uniform_int_distribution<int> coef(-6, 6);
  const int nb = ni + nc;
  std::vector<double> cost(nb), lb(nb), ub(nb);
  for (int j = 0; j < nb; ++j) {
    cost[j] = coef(rng);
    lb[j] = 0;
    ub[j] = j < ni ? 1 + rng() % 3 : 2.5;
  }
  struct Row {
    double lo, up;
    std::vector<std::pair<int, double>> t;  // (block column or -1 - global column, value)
  };
  std::vector<Row> block(nr);
  for (auto& r : block) {
    double act = 0.0;
    for (int j = 0; j < nb; ++j) {
      if (rng() % 2) continue;
      const double a = coef(rng) + (j >= ni ? 0.5 : 0.0);
      if (a == 0.0) continue;
      r.t.emplace_back(j, a);
      act += a * std::floor(ub[j] / 2);
    }
    r.lo = rng() % 2 ? -kInf : act - 2;
    r.up = act + 1 + rng() % 3;
  }
  // Global column order: a random permutation of all (copy, block column).
  std::vector<std::pair<int, int>> cols;
  for (int c = 0; c < k; ++c) {
    for (int j = 0; j < nb; ++j) cols.emplace_back(c, j);
  }
  std::shuffle(cols.begin(), cols.end(), rng);
  std::map<std::pair<int, int>, Index> id;
  ModelBuilder b;
  for (auto [c, j] : cols) {
    id[{c, j}] = b.add_col("c" + std::to_string(c) + "_" + std::to_string(j), cost[j], lb[j], ub[j],
                           j < ni ? VarType::Integer : VarType::Continuous);
  }
  std::vector<std::pair<int, int>> rows;  // (copy or -1 for linking, row)
  for (int c = 0; c < k; ++c) {
    for (int r = 0; r < nr; ++r) rows.emplace_back(c, r);
  }
  for (int r = 0; r < nlink; ++r) rows.emplace_back(-1, r);
  std::shuffle(rows.begin(), rows.end(), rng);
  std::vector<std::vector<std::pair<int, double>>> link(nlink);
  std::vector<double> link_up(nlink);
  for (int r = 0; r < nlink; ++r) {
    for (int j = 0; j < nb; ++j) {
      if (rng() % 3 == 0) link[r].emplace_back(j, 1 + rng() % 3);
    }
    link_up[r] = static_cast<double>(k * (1 + rng() % 2));
  }
  int cnt = 0;
  for (auto [c, r] : rows) {
    if (c >= 0) {
      const Index row = b.add_row("r" + std::to_string(cnt++), block[r].lo, block[r].up);
      for (auto [j, a] : block[r].t) b.add_coef(row, id[{c, j}], a);
    } else {
      const Index row = b.add_row("l" + std::to_string(cnt++), -kInf, link_up[r]);
      for (int cc = 0; cc < k; ++cc) {
        for (auto [j, a] : link[r]) b.add_coef(row, id[{cc, j}], a);
      }
    }
  }
  return b.build();
}

MipProblem to_problem(const Model& m) {
  MipProblem p;
  MipProblem::from_model(m, &p);
  return p;
}

double factorial(int k) {
  double f = 1.0;
  for (int i = 2; i <= k; ++i) f *= i;
  return f;
}

}  // namespace

TEST(symmetry_small_graph_groups) {
  // Petersen graph: |Aut| = 120, vertex transitive.
  std::vector<std::pair<int, int>> pet = cycle(5);
  for (auto [a, b] : std::vector<std::pair<int, int>>{{5, 7}, {7, 9}, {9, 6}, {6, 8}, {8, 5}}) pet.emplace_back(a, b);
  for (int i = 0; i < 5; ++i) pet.emplace_back(i, i + 5);
  {
    const MipProblem p = to_problem(independent_set(10, pet));
    const Symmetry s = detect_symmetry(p);
    CHECK(all_generators_valid(p, s));
    CHECK(!s.aborted);
    CHECK_NEAR(group_order(s), 120.0, 0.5);
    CHECK(full_orbit(s, 0, 10).size() == 10u);
  }
  // 3-cube: |Aut| = 48.
  std::vector<std::pair<int, int>> cube;
  for (int v = 0; v < 8; ++v) {
    for (int b = 0; b < 3; ++b) {
      if (v < (v ^ (1 << b))) cube.emplace_back(v, v ^ (1 << b));
    }
  }
  {
    const MipProblem p = to_problem(independent_set(8, cube));
    const Symmetry s = detect_symmetry(p);
    CHECK(all_generators_valid(p, s));
    CHECK_NEAR(group_order(s), 48.0, 0.5);
  }
  // C6 plus two triangles: all vertices are 2-regular, so refinement alone
  // cannot separate them; the group is D6 x (S3 wr S2) of order 12 * 72.
  {
    std::vector<std::pair<int, int>> e = cycle(6);
    for (auto x : cycle(3, 6)) e.push_back(x);
    for (auto x : cycle(3, 9)) e.push_back(x);
    const MipProblem p = to_problem(independent_set(12, e));
    const Symmetry s = detect_symmetry(p);
    CHECK(all_generators_valid(p, s));
    CHECK_NEAR(group_order(s), 864.0, 0.5);
    CHECK(full_orbit(s, 0, 12) == (std::vector<Index>{0, 1, 2, 3, 4, 5}));
    CHECK(full_orbit(s, 7, 12) == (std::vector<Index>{6, 7, 8, 9, 10, 11}));
  }
  // Weights break the symmetry of a cycle down to a reflection.
  {
    const MipProblem p = to_problem(independent_set(6, cycle(6), {1, 2, 3, 4, 3, 2}));
    const Symmetry s = detect_symmetry(p);
    CHECK(all_generators_valid(p, s));
    CHECK_NEAR(group_order(s), 2.0, 0.1);
    CHECK(full_orbit(s, 0, 6).size() == 1u);
    CHECK(full_orbit(s, 1, 6) == (std::vector<Index>{1, 5}));
  }
}

TEST(symmetry_bin_packing_and_coloring) {
  {
    // 4 identical bins x 2 equal items: S4 x S2.
    const MipProblem p = to_problem(bin_packing({3, 3, 4, 5, 6}, 4, 10));
    const Symmetry s = detect_symmetry(p);
    CHECK(all_generators_valid(p, s));
    CHECK_NEAR(group_order(s), 48.0, 0.5);
    CHECK(full_orbit(s, 0, p.num_cols()) == (std::vector<Index>{0, 1, 2, 3}));
  }
  {
    // Coloring a 5-cycle with 4 colors: D5 x S4.
    const MipProblem p = to_problem(coloring(5, cycle(5), 4));
    const Symmetry s = detect_symmetry(p);
    CHECK(all_generators_valid(p, s));
    CHECK_NEAR(group_order(s), 10.0 * 24.0, 0.5);
  }
}

TEST(symmetry_none_when_asymmetric) {
  std::mt19937 rng(3);
  for (int t = 0; t < 20; ++t) {
    ModelBuilder b;
    const int n = 6 + t;
    for (int j = 0; j < n; ++j) b.add_col("x" + std::to_string(j), j + 1.0, 0, 1, VarType::Binary);
    for (int r = 0; r < 4; ++r) {
      const Index row = b.add_row("r" + std::to_string(r), -kInf, 3);
      for (int j = 0; j < n; ++j) {
        if (rng() % 2) b.add_coef(row, j, 1 + rng() % 3);
      }
    }
    const MipProblem p = to_problem(b.build());
    const Symmetry s = detect_symmetry(p);
    CHECK(s.empty());
  }
}

TEST(symmetry_stabilizer_chain_valid) {
  // Repeated stabilizers (as along a path of up-branchings) stay automorphisms
  // and fix every point stabilized so far; on S_k-like groups they stay large.
  std::mt19937 rng(7);
  GroupScratch scratch;
  for (int t = 0; t < scaled(30); ++t) {
    Model m;
    if (t % 3 == 0) {
      m = coloring(6 + t % 4, cycle(6 + t % 4), 5);
    } else if (t % 3 == 1) {
      std::vector<std::pair<int, int>> pet = cycle(5);
      for (auto [a, b] : std::vector<std::pair<int, int>>{{5, 7}, {7, 9}, {9, 6}, {6, 8}, {8, 5}}) pet.emplace_back(a, b);
      for (int i = 0; i < 5; ++i) pet.emplace_back(i, i + 5);
      m = independent_set(10, pet);
    } else {
      m = duplicated_blocks(rng, 3 + rng() % 3, 2 + rng() % 3, 0, 1 + rng() % 2, 1);
    }
    const MipProblem p = to_problem(m);
    const Index n = p.num_cols();
    const Symmetry s = detect_symmetry(p);
    PermGroup g;
    g.gens = s.gens;
    std::vector<Index> fixed;
    std::vector<Index> orb;
    std::vector<char> mark(n, 0);
    for (int step = 0; step < 6 && !g.gens.empty(); ++step) {
      // Branch on a moved column, as the search would.
      Index j = -1;
      for (int tries = 0; tries < 50 && j < 0; ++tries) {
        const Index c = static_cast<Index>(rng() % n);
        g.orbit(c, orb, mark);
        if (orb.size() > 1) j = c;
      }
      if (j < 0) break;
      g = g.stabilizer(j, n, rng(), scratch);
      fixed.push_back(j);
      for (size_t k = 0; k < g.gens.size(); ++k) {
        const std::vector<Index> sigma = [&] {
          std::vector<Index> v(n);
          for (Index c = 0; c < n; ++c) v[c] = c;
          for (const auto& [a, b] : g.gens[k]) v[a] = b;
          return v;
        }();
        CHECK(maps_onto_itself(p, sigma));
        for (Index f : fixed) CHECK(sigma[f] == f);
        CHECK(std::is_sorted(g.gens[k].begin(), g.gens[k].end()));
      }
    }
  }
  {
    // Coloring with 6 identical colors: after fixing x[0][c], colors other
    // than c stay interchangeable (orbit of w_c' has 5 elements).
    const MipProblem p = to_problem(coloring(4, cycle(4), 6));
    const Symmetry s = detect_symmetry(p);
    PermGroup g;
    g.gens = s.gens;
    const Index x00 = 6;  // columns: w_0..w_5, then x[v][c]
    g = g.stabilizer(x00, p.num_cols(), 1, scratch);
    std::vector<Index> orb;
    std::vector<char> mark(p.num_cols(), 0);
    g.orbit(1, orb, mark);
    CHECK(orb.size() == 5u);
    g.orbit(0, orb, mark);
    CHECK(orb.size() == 1u);
  }
}

TEST(symmetry_duplicated_blocks_detected) {
  std::mt19937 rng(11);
  for (int t = 0; t < scaled(100); ++t) {
    const int k = 2 + rng() % 4;
    const MipProblem p = to_problem(duplicated_blocks(rng, k, 2 + rng() % 4, rng() % 2, 1 + rng() % 3, rng() % 3));
    const Symmetry s = detect_symmetry(p);
    CHECK(all_generators_valid(p, s));
    // The copies are interchangeable: the group contains S_k.
    CHECK(group_order(s) >= factorial(k) - 0.5);
    // Each generator really moves something and the chain is consistent.
    for (int g = 0; g < s.num_generators(); ++g) {
      const std::vector<Index> sigma = generator(s, g, p.num_cols());
      for (int l = 0; l < s.gen_level[g]; ++l) CHECK(sigma[s.base[l]] == s.base[l]);
    }
    for (size_t l = 0; l < s.tree.size(); ++l) {
      for (const auto& e : s.tree[l]) {
        std::vector<Index> tr;
        REQUIRE(s.transversal(static_cast<int>(l), e.point, tr));
        CHECK(tr[s.base[l]] == e.point);
        CHECK(maps_onto_itself(p, tr));
      }
    }
  }
}

namespace {

struct Outcome {
  Status status;
  double obj, bound;
  bool feasible;
  long long nodes;
};

Outcome run(const Model& m, bool sym, int threads, bool presolve, bool det) {
  set_env("MILP_NO_SYMMETRY", !sym);
  SolverParams sp;
  sp.verbosity = 0;
  sp.threads = threads;
  sp.presolve = presolve;
  sp.deterministic = det;
  sp.time_limit = 2;
  sp.mip_rel_gap = 0.0;
  sp.mip_abs_gap = 1e-9;
  const SolveResult r = solve(m, sp);
  set_env("MILP_NO_SYMMETRY", false);
  bool feas = true;
  if (r.has_solution) feas = m.violation(r.x).max() <= 1e-5;
  return {r.status, r.has_solution ? r.objective : kInf, r.bound, feas, r.stats.nodes};
}

// Same optimum with and without symmetry handling; if one side times out,
// its incumbent and bound must still be consistent with the other's optimum.
void compare(const Model& m, int variant, int& checked, long long& branched) {
  const int threads = variant % 3 == 0 ? 1 : variant % 3 == 1 ? 2 : 3;
  const bool presolve = variant % 2 == 0, det = variant % 3 == 2;
  const Outcome a = run(m, false, threads, presolve, det);
  const Outcome b = run(m, true, threads, presolve, det);
  CHECK(b.feasible);
  // Orbital branching alone (orbitopal fixing and symmetry-breaking rows off) on the same model.
  set_env("MILP_NO_ORBITOPE", true);
  set_env_value("MILP_SST", "0");
  const Outcome c = run(m, true, threads, presolve, det);
  set_env("MILP_NO_ORBITOPE", false);
  set_env("MILP_SST", false);
  CHECK(c.feasible);
  branched += c.nodes > 1;
  if (a.status == Status::Optimal && c.status == Status::Optimal) {
    CHECK(std::abs(a.obj - c.obj) <= 1e-6 * (1.0 + std::abs(a.obj)));
  }
  if (std::getenv("MILP_TEST_VERBOSE")) {
    std::fprintf(stderr, "  %s v%d: obj %g / %g, bound %g / %g, nodes %lld / %lld\n", m.name.c_str(), variant, a.obj,
                 b.obj, a.bound, b.bound, a.nodes, b.nodes);
  }
  const double tol = 1e-6 * (1.0 + std::abs(a.obj == kInf ? b.obj : a.obj));
  if (a.status == Status::Optimal && b.status == Status::Optimal) {
    ++checked;
    if (std::abs(a.obj - b.obj) > tol) {
      std::fprintf(stderr, "  %s: optimum %.10g without symmetry, %.10g with (variant %d)\n", m.name.c_str(), a.obj,
                   b.obj, variant);
      CHECK(false);
    }
  } else if (a.status == Status::Optimal && b.status == Status::TimeLimit) {
    CHECK(b.obj >= a.obj - tol);
    CHECK(b.bound <= a.obj + tol);
  } else if (b.status == Status::Optimal && a.status == Status::TimeLimit) {
    CHECK(a.obj >= b.obj - tol);
    CHECK(a.bound <= b.obj + tol);
  } else if (a.status == Status::Infeasible || b.status == Status::Infeasible) {
    CHECK(a.status == b.status);
  }
}

}  // namespace

TEST(symmetry_end_to_end_same_optimum) {
  std::mt19937 rng(2024);
  int checked = 0;
  long long branched = 0;
  // Zero-half cuts (odd cycles) solve most independent set instances at the root; the
  // test needs symmetric searches that branch.
  set_env("MILP_NO_ZEROHALF", true);
  for (int t = 0; t < scaled(30); ++t) {
    Model m = duplicated_blocks(rng, 2 + rng() % 4, 2 + rng() % 4, rng() % 2, 1 + rng() % 3, 1 + rng() % 2);
    m.name = "blocks" + std::to_string(t);
    compare(m, t, checked, branched);
  }
  for (int t = 0; t < scaled(100); ++t) {
    std::vector<int> sizes;
    const int n = 8 + rng() % 5;
    for (int i = 0; i < n; ++i) sizes.push_back(4 + rng() % 9);
    Model m = bin_packing(sizes, 2 + n / 2, 17 + rng() % 4);
    m.name = "binpack" + std::to_string(t);
    compare(m, t, checked, branched);
  }
  for (int t = 0; t < scaled(40); ++t) {
    // Coloring of random graphs (identical colors) and independent sets on
    // circulant graphs (vertex transitive).
    std::vector<std::pair<int, int>> e;
    if (t % 2 == 0) {
      const int n = 10 + rng() % 6;
      for (int u = 0; u < n; ++u) {
        for (int v = u + 1; v < n; ++v) {
          if (rng() % 2) e.emplace_back(u, v);
        }
      }
      Model m = coloring(n, e, 5 + rng() % 2);
      m.name = "coloring" + std::to_string(t);
      compare(m, t, checked, branched);
    } else {
      const int n = 14 + rng() % 12;
      for (int d = 1; d <= n / 2; ++d) {
        if (d == 1 || rng() % 3 == 0) {
          for (int i = 0; i < n; ++i) {
            if (d * 2 != n || i < d) e.emplace_back(i, (i + d) % n);
          }
        }
      }
      Model m = independent_set(n, e);
      m.name = "indset" + std::to_string(t);
      compare(m, t, checked, branched);
    }
  }
  set_env("MILP_NO_ZEROHALF", false);
  CHECK(checked > 30);
  // The symmetric searches really branch (12-15 of the instances; the count moves with LP pivoting
  // details, so the bound leaves some room).
  CHECK(branched >= 10);
}

TEST(orbitope_detected_on_coloring_and_bin_packing) {
  {
    // Coloring a random graph with 5 identical colors: rows "vertex v gets one color" are packing rows.
    std::vector<std::pair<int, int>> e{{0, 1}, {1, 2}, {2, 3}, {3, 0}, {0, 2}, {4, 5}, {5, 6}, {6, 4}, {1, 5}};
    const MipProblem p = to_problem(coloring(7, e, 5));
    const Symmetry s = detect_symmetry(p);
    std::vector<char> used;
    const std::vector<Orbitope> orbs = detect_packing_orbitopes(p, s, &used);
    CHECK(orbs.size() == 1);
    if (!orbs.empty()) {
      CHECK(orbs[0].cols == 5);
      CHECK(orbs[0].rows >= 1);
      // Swapping two matrix columns (colors with all their columns) is a symmetry.
      // (Checked by the detection itself; here: entries are distinct binaries.)
      std::vector<Index> v = orbs[0].var;
      std::sort(v.begin(), v.end());
      CHECK(std::unique(v.begin(), v.end()) == v.end());
    }
  }
  {
    const MipProblem p = to_problem(bin_packing({3, 4, 5, 6, 7}, 4, 10));
    const Symmetry s = detect_symmetry(p);
    const std::vector<Orbitope> orbs = detect_packing_orbitopes(p, s, nullptr);
    for (const Orbitope& o : orbs) CHECK(o.cols >= 2 && o.rows >= 1);
  }
}

TEST(orbitope_fixing_keeps_all_representatives) {
  // Enumerate all r x p packing matrices (each row at most one 1); the representatives have strictly
  // increasing leading rows of the nonzero columns, zero columns last. Under random fixings, the
  // propagation must keep every representative that satisfies them and fail only if none does.
  std::mt19937 rng(7);
  for (int t = 0; t < scaled(400); ++t) {
    const Index r = 1 + rng() % 4, pc = 1 + rng() % 4;
    Orbitope o;
    o.rows = r, o.cols = pc;
    for (Index k = 0; k < r * pc; ++k) o.var.push_back(k);
    std::vector<double> lb(r * pc, 0.0), ub(r * pc, 1.0);
    for (Index k = 0; k < r * pc; ++k) {
      const int f = rng() % 6;
      if (f == 0) ub[k] = 0.0;
      if (f == 1) lb[k] = 1.0;
    }
    std::vector<double> lb0 = lb, ub0 = ub;
    std::vector<Index> tight;
    std::vector<char> ones;
    const bool ok = orbitope_propagate(o, lb, ub, tight, ones);
    // Enumerate: each row picks a column or none ((pc + 1)^r).
    long long total = 1;
    for (Index i = 0; i < r; ++i) total *= pc + 1;
    bool any = false;
    for (long long code = 0; code < total; ++code) {
      std::vector<int> x(r * pc, 0);
      long long c = code;
      for (Index i = 0; i < r; ++i) {
        const int col = static_cast<int>(c % (pc + 1));
        c /= pc + 1;
        if (col < pc) x[i * pc + col] = 1;
      }
      bool fits = true;
      for (Index k = 0; k < r * pc; ++k) fits = fits && x[k] >= lb0[k] && x[k] <= ub0[k];
      if (!fits) continue;
      // Representative check.
      int prev = -1;
      bool rep = true, zero_seen = false;
      for (Index col = 0; col < pc; ++col) {
        int lead = -1;
        for (Index i = 0; i < r && lead < 0; ++i) {
          if (x[i * pc + col]) lead = i;
        }
        if (lead < 0) {
          zero_seen = true;
          continue;
        }
        if (zero_seen || lead <= prev) rep = false;
        prev = lead;
      }
      if (!rep) continue;
      any = true;
      CHECK(ok);
      for (Index k = 0; k < r * pc && ok; ++k) CHECK(x[k] >= lb[k] && x[k] <= ub[k]);
    }
    if (!any) CHECK(!ok || true);  // failing early is optional, keeping a representative is not
  }
}

TEST(orbitope_canonical_solutions_survive_fixing) {
  // Optimal colorings mapped onto their representative satisfy orbitopal fixing at full assignment.
  std::mt19937 rng(11);
  int tested = 0;
  for (int t = 0; t < scaled(10); ++t) {
    std::vector<std::pair<int, int>> e;
    const int nv = 6 + rng() % 4;
    for (int u = 0; u < nv; ++u)
      for (int v = u + 1; v < nv; ++v)
        if (rng() % 3 == 0) e.emplace_back(u, v);
    const Model m = coloring(nv, e, 4);
    const MipProblem p = to_problem(m);
    const Symmetry s = detect_symmetry(p);
    const std::vector<Orbitope> orbs = detect_packing_orbitopes(p, s, nullptr);
    SolverParams prm;
    prm.verbosity = 0;
    prm.threads = 1;
    prm.presolve = false;
    const SolveResult r = solve(m, prm);
    if (r.status != Status::Optimal || orbs.empty()) continue;
    std::vector<double> x = r.x;
    for (const Orbitope& o : orbs) orbitope_canonicalize(o, x);
    CHECK(p.max_violation(x) <= 1e-6);
    CHECK(std::abs(p.objective(x) - p.objective(r.x)) <= 1e-6);
    for (const Orbitope& o : orbs) {
      std::vector<double> lb = x, ub = x;
      std::vector<Index> tight;
      std::vector<char> ones;
      CHECK(orbitope_propagate(o, lb, ub, tight, ones));
      CHECK(tight.empty());
    }
    ++tested;
  }
  CHECK(tested >= 3);
}

namespace {

// Number of column permutations mapping p onto itself (brute force, n <= 8).
long long count_automorphisms(const MipProblem& p) {
  const Index n = p.num_cols();
  std::vector<Index> sigma(n);
  for (Index j = 0; j < n; ++j) sigma[j] = j;
  long long cnt = 0;
  do {
    cnt += maps_onto_itself(p, sigma);
  } while (std::next_permutation(sigma.begin(), sigma.end()));
  return cnt;
}

// Optimum of a pure binary problem by enumeration (kInf if infeasible).
double enumerate_binary(const MipProblem& p) {
  const Index n = p.num_cols(), m = p.num_rows();
  const SparseMatrix R = p.A.transpose();
  double best = kInf;
  std::vector<double> x(n);
  for (long long mask = 0; mask < (1LL << n); ++mask) {
    bool ok = true;
    for (Index j = 0; j < n && ok; ++j) {
      x[j] = static_cast<double>((mask >> j) & 1);
      ok = x[j] >= p.col_lb[j] && x[j] <= p.col_ub[j];
    }
    for (Index r = 0; r < m && ok; ++r) {
      double a = 0.0;
      for (Index k = R.begin(r); k < R.end(r); ++k) a += R.value[k] * x[R.index[k]];
      ok = a >= p.row_lb[r] - 1e-9 && a <= p.row_ub[r] + 1e-9;
    }
    if (ok) best = std::min(best, p.objective(x));
  }
  return best;
}

// Small pure binary models with symmetry: independent sets on vertex-transitive graphs (weights
// equal or equal on orbits of a rotation), colorings, bin packing, a covering design.
Model small_symmetric_binary(std::mt19937& rng, int t) {
  switch (t % 5) {
    case 0: {
      const int n = 8 + rng() % 7;
      std::vector<std::pair<int, int>> e;
      for (int d = 1; d <= n / 2; ++d) {
        if (d == 1 || rng() % 3 == 0) {
          for (int i = 0; i < n; ++i) {
            if (d * 2 != n || i < d) e.emplace_back(i, (i + d) % n);
          }
        }
      }
      std::vector<double> w;
      if (rng() % 2) {
        // Weights invariant under the rotation by 2 (n even) or none.
        for (int i = 0; i < n; ++i) w.push_back(n % 2 == 0 ? 1.0 + (i % 2) : 1.0);
      }
      return independent_set(n, e, w);
    }
    case 1: {
      // Hypercube Q4 with its antipodal pairs also in conflict, or the Petersen graph.
      if (rng() % 2) {
        std::vector<std::pair<int, int>> e;
        for (int v = 0; v < 16; ++v) {
          for (int b = 0; b < 4; ++b) {
            if (v < (v ^ (1 << b))) e.emplace_back(v, v ^ (1 << b));
          }
          if (v < 15 - v) e.emplace_back(v, 15 - v);
        }
        return independent_set(16, e);
      }
      std::vector<std::pair<int, int>> pet = cycle(5);
      for (auto [a, b] : std::vector<std::pair<int, int>>{{5, 7}, {7, 9}, {9, 6}, {6, 8}, {8, 5}}) pet.emplace_back(a, b);
      for (int i = 0; i < 5; ++i) pet.emplace_back(i, i + 5);
      return independent_set(10, pet);
    }
    case 2: {
      std::vector<std::pair<int, int>> e;
      const int nv = 3 + rng() % 2;
      for (int u = 0; u < nv; ++u)
        for (int v = u + 1; v < nv; ++v)
          if (rng() % 3) e.emplace_back(u, v);
      return coloring(nv, e, 3);
    }
    case 3: {
      std::vector<int> sizes;
      for (int i = 0; i < 4; ++i) sizes.push_back(2 + rng() % 3);
      return bin_packing(sizes, 3, 6);
    }
    default: {
      // Circulant covering (each triple i, i + d, i + 2d + 1 has a chosen column) with a
      // cardinality row: symmetric under the rotations.
      const int n = 9 + rng() % 6;
      ModelBuilder b;
      for (int j = 0; j < n; ++j) b.add_col("x" + std::to_string(j), 1, 0, 1, VarType::Binary);
      const int d = 1 + rng() % 3;
      for (int i = 0; i < n; ++i) {
        const Index r = b.add_row("c" + std::to_string(i), 1, kInf);
        b.add_coef(r, i, 1);
        b.add_coef(r, (i + d) % n, 1);
        b.add_coef(r, (i + 2 * d + 1) % n, 1);
      }
      const Index r = b.add_row("card", 2, n - 1);
      for (int j = 0; j < n; ++j) b.add_coef(r, j, 1);
      return b.build();
    }
  }
}

}  // namespace

TEST(symmetry_group_order_matches_brute_force) {
  // The automorphism search on small random symmetric models: the group order from the
  // stabilizer chain equals the number of column automorphisms found by brute force.
  std::mt19937 rng(5);
  int compared = 0;
  for (int t = 0; t < scaled(60); ++t) {
    ModelBuilder b;
    const int n = 3 + rng() % 6;
    std::vector<double> cost(n);
    for (int j = 0; j < n; ++j) cost[j] = static_cast<double>(rng() % 2);
    for (int j = 0; j < n; ++j) b.add_col("x" + std::to_string(j), cost[j], 0, 1, VarType::Binary);
    // Random graph-like rows (pairs with coefficient 1 or 2), so that there is often symmetry.
    const int nr = 1 + rng() % (2 * n);
    for (int r = 0; r < nr; ++r) {
      const Index row = b.add_row("r" + std::to_string(r), -kInf, 1 + rng() % 2);
      const int u = rng() % n, v = (u + 1 + rng() % (n - 1)) % n;
      const double c = 1.0 + rng() % 2;
      b.add_coef(row, u, c);
      b.add_coef(row, v, c);
    }
    const MipProblem p = to_problem(b.build());
    const Symmetry s = detect_symmetry(p);
    CHECK(all_generators_valid(p, s));
    if (s.aborted) continue;
    const long long bf = count_automorphisms(p);
    CHECK_NEAR(group_order(s), static_cast<double>(bf), 0.5);
    ++compared;
  }
  CHECK(compared >= scaled(50));
  // Complete graph K6: S6.
  std::vector<std::pair<int, int>> k6;
  for (int u = 0; u < 6; ++u)
    for (int v = u + 1; v < 6; ++v) k6.emplace_back(u, v);
  const MipProblem p = to_problem(independent_set(6, k6));
  CHECK_NEAR(group_order(detect_symmetry(p)), 720.0, 0.5);
}

TEST(sst_rows_keep_the_optimum) {
  // Symmetry-breaking rows and fixings (search/sst.hpp) on small symmetric binary models: the
  // optimum by enumeration is unchanged, for both modes, several seeds and with the generators of
  // packing orbitopes skipped or not.
  std::mt19937 rng(17);
  int reduced = 0, fixed = 0;
  for (int t = 0; t < scaled(60); ++t) {
    const MipProblem p = to_problem(small_symmetric_binary(rng, t));
    if (p.num_cols() > 18) continue;
    const Symmetry s = detect_symmetry(p);
    CHECK(all_generators_valid(p, s));
    const double opt = enumerate_binary(p);
    std::vector<char> used;
    detect_packing_orbitopes(p, s, &used);
    for (int variant = 0; variant < 4; ++variant) {
      SstOptions o;
      o.mode = 1 + variant % 2;
      o.seed = static_cast<uint64_t>(t * 7 + variant);
      const SstResult r = sst_reductions(p, s, variant < 2 ? std::vector<char>{} : used, o);
      // Rows only between columns of one orbit (equal data), fixings only of binaries.
      for (const auto& [l, f] : r.rows) {
        CHECK(l != f);
        CHECK(p.cost[l] == p.cost[f] && p.col_lb[l] == p.col_lb[f] && p.col_ub[l] == p.col_ub[f]);
      }
      if (r.empty()) continue;
      ++reduced;
      fixed += !r.fixed_zero.empty();
      const MipProblem q = apply_sst(p, r);
      const double opt2 = enumerate_binary(q);
      if (std::abs(opt2 - opt) > 1e-9) {
        std::fprintf(stderr, "  model %d variant %d: optimum %g, with SST rows %g (%zu rows, %zu fixed)\n", t, variant,
                     opt, opt2, r.rows.size(), r.fixed_zero.size());
        CHECK(false);
      }
    }
  }
  CHECK(reduced >= scaled(60));
  CHECK(fixed >= 10);
}

TEST(symmetry_solves_match_enumeration) {
  // End to end with all symmetry handling (SST rows, orbitopal fixing, orbital branching) against
  // enumeration, opportunistic with 1 / 3 threads and deterministic with 2.
  std::mt19937 rng(23);
  set_env("MILP_NO_ZEROHALF", true);
  set_env_value("MILP_SST_MIN_FIX", "0");  // symmetry-breaking rows whenever they fix something
  int checked = 0;
  for (int t = 0; t < scaled(40); ++t) {
    const Model m = small_symmetric_binary(rng, t);
    const MipProblem p = to_problem(m);
    if (p.num_cols() > 18) continue;
    const double opt = enumerate_binary(p);
    for (int variant = 0; variant < 3; ++variant) {
      SolverParams sp;
      sp.verbosity = 0;
      sp.threads = variant == 0 ? 1 : variant == 1 ? 3 : 2;
      sp.deterministic = variant == 2;
      sp.time_limit = 5;
      sp.mip_rel_gap = 0.0;
      sp.mip_abs_gap = 1e-9;
      const SolveResult r = solve(m, sp);
      if (opt == kInf) {
        CHECK(r.status == Status::Infeasible);
        continue;
      }
      REQUIRE(r.status == Status::Optimal);
      CHECK(m.violation(r.x).max() <= 1e-6);
      if (std::abs(r.objective - opt) > 1e-6) {
        std::fprintf(stderr, "  model %d variant %d: optimum %g, solve %g\n", t, variant, opt, r.objective);
        CHECK(false);
      }
      ++checked;
    }
  }
  set_env("MILP_NO_ZEROHALF", false);
  set_env("MILP_SST_MIN_FIX", false);
  CHECK(checked >= scaled(90));
}
