// Tests for the anonymized model export (io/anonymize.hpp, popso --write-anon) and the MPS
// writer's exact numbers: anonymize, write, read back, compare the data exactly (as multisets and
// through the map) and solve both.

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "milp/io/anonymize.hpp"
#include "milp/io/mps_reader.hpp"
#include "milp/solver/solver.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

// Exact equality (bit for bit, except that 0.0 and -0.0 are the same: the writer leaves zeros out).
bool same(double a, double b) { return a == b; }

using Entries = std::vector<std::pair<Index, double>>;  // (column, value) of a row

std::vector<Entries> rows_of(const Model& m) {
  std::vector<Entries> rows(static_cast<size_t>(m.num_rows()));
  for (Index j = 0; j < m.num_cols(); ++j)
    for (Index p = m.A.begin(j); p < m.A.end(j); ++p)
      if (m.A.value[p] != 0.0) rows[m.A.index[p]].emplace_back(j, m.A.value[p]);
  return rows;
}

// Column key (obj, lb, ub, type) and row key (lb, ub, sorted (coefficient, column key) pairs):
// permutation-free descriptions whose multisets must agree between a model and its copy.
using Key = std::vector<double>;

Key col_key(const Model& m, Index j) {
  return {m.obj[j], m.col_lb[j], m.col_ub[j], static_cast<double>(static_cast<int>(m.col_type[j]))};
}

std::vector<Key> col_multiset(const Model& m) {
  std::vector<Key> keys;
  for (Index j = 0; j < m.num_cols(); ++j) keys.push_back(col_key(m, j));
  std::sort(keys.begin(), keys.end());
  return keys;
}

// Rows with a finite side (free rows are written as extra N rows, which readers drop).
std::vector<Key> row_multiset(const Model& m) {
  std::vector<Key> keys;
  const std::vector<Entries> rows = rows_of(m);
  for (Index i = 0; i < m.num_rows(); ++i) {
    if (m.row_lb[i] == -kInf && m.row_ub[i] == kInf) continue;
    std::vector<Key> entries;
    for (const auto& [j, v] : rows[i]) {
      Key e = col_key(m, j);
      e.insert(e.begin(), v);
      entries.push_back(e);
    }
    std::sort(entries.begin(), entries.end());
    Key k = {m.row_lb[i], m.row_ub[i]};
    for (const Key& e : entries) k.insert(k.end(), e.begin(), e.end());
    keys.push_back(k);
  }
  std::sort(keys.begin(), keys.end());
  return keys;
}

bool is_permutation_of(const std::vector<Index>& p, size_t n) {
  if (p.size() != n) return false;
  std::vector<char> seen(n, 0);
  for (Index k : p) {
    if (k < 0 || static_cast<size_t>(k) >= n || seen[static_cast<size_t>(k)]) return false;
    seen[static_cast<size_t>(k)] = 1;
  }
  return true;
}

// The copy `a` (as anonymize returned it, or as read back from its file) is `orig` permuted by
// `map`: every column, row, SOS set and indicator constraint, entry by entry, exactly.
void check_mapped(const Model& orig, const Model& a, const AnonMap& map, bool free_rows_dropped) {
  REQUIRE(is_permutation_of(map.col, static_cast<size_t>(orig.num_cols())));
  REQUIRE(is_permutation_of(map.row, static_cast<size_t>(orig.num_rows())));
  REQUIRE(is_permutation_of(map.sos, orig.sos.size()));
  REQUIRE(is_permutation_of(map.indicator, orig.indicators.size()));
  REQUIRE(a.num_cols() == orig.num_cols());
  CHECK(a.sense == orig.sense);
  CHECK(same(a.obj_offset, orig.obj_offset));
  std::vector<Index> new_col(map.col.size());
  for (size_t k = 0; k < map.col.size(); ++k) new_col[static_cast<size_t>(map.col[k])] = static_cast<Index>(k);
  for (Index j = 0; j < a.num_cols(); ++j) CHECK(col_key(a, j) == col_key(orig, map.col[j]));
  // Rows: the copy's rows in order are the original's rows map.row (free ones left out if dropped).
  std::vector<Index> expect_rows;
  for (Index o : map.row)
    if (!free_rows_dropped || orig.row_lb[o] > -kInf || orig.row_ub[o] < kInf) expect_rows.push_back(o);
  REQUIRE(a.num_rows() == static_cast<Index>(expect_rows.size()));
  const std::vector<Entries> ra = rows_of(a), ro = rows_of(orig);
  for (Index i = 0; i < a.num_rows(); ++i) {
    const Index o = expect_rows[i];
    CHECK(same(a.row_lb[i], orig.row_lb[o]));
    CHECK(same(a.row_ub[i], orig.row_ub[o]));
    Entries mapped;
    for (const auto& [j, v] : ro[o]) mapped.emplace_back(new_col[j], v);
    std::sort(mapped.begin(), mapped.end());
    Entries got = ra[i];
    std::sort(got.begin(), got.end());
    CHECK(got == mapped);
  }
  // Quadratic objective: the assembled lower triangles agree under the permutation.
  {
    const SparseMatrix qa = a.q_lower(), qo = orig.q_lower();
    std::vector<std::tuple<Index, Index, double>> ea, eo;
    for (Index j = 0; j < qa.num_major; ++j)
      for (Index p = qa.begin(j); p < qa.end(j); ++p) ea.emplace_back(qa.index[p], j, qa.value[p]);
    for (Index j = 0; j < qo.num_major; ++j) {
      for (Index p = qo.begin(j); p < qo.end(j); ++p) {
        const Index r = new_col[qo.index[p]], c = new_col[j];
        eo.emplace_back(std::max(r, c), std::min(r, c), qo.value[p]);
      }
    }
    std::sort(ea.begin(), ea.end());
    std::sort(eo.begin(), eo.end());
    CHECK(ea == eo);
    for (size_t k = 0; k < a.q_row.size(); ++k) CHECK(a.q_row[k] >= a.q_col[k]);
  }
  REQUIRE(a.sos.size() == orig.sos.size());
  for (size_t s = 0; s < a.sos.size(); ++s) {
    const SosConstraint& x = a.sos[s];
    const SosConstraint& y = orig.sos[static_cast<size_t>(map.sos[s])];
    CHECK(x.type == y.type);
    REQUIRE(x.cols.size() == y.cols.size());
    for (size_t k = 0; k < x.cols.size(); ++k) {
      CHECK(x.cols[k] == new_col[y.cols[k]]);
      CHECK(same(x.weights[k], y.weights[k]));
    }
  }
  // Indicator constraints (all with a finite side here: the writer leaves vacuous ones out).
  REQUIRE(a.indicators.size() == orig.indicators.size());
  for (size_t q = 0; q < a.indicators.size(); ++q) {
    const IndicatorConstraint& x = a.indicators[q];
    const IndicatorConstraint& y = orig.indicators[static_cast<size_t>(map.indicator[q])];
    CHECK(x.binary == new_col[y.binary]);
    CHECK(x.value == y.value);
    CHECK(same(x.lb, y.lb));
    CHECK(same(x.ub, y.ub));
    Entries ex, ey;
    for (size_t k = 0; k < x.cols.size(); ++k) ex.emplace_back(x.cols[k], x.coefs[k]);
    for (size_t k = 0; k < y.cols.size(); ++k) ey.emplace_back(new_col[y.cols[k]], y.coefs[k]);
    std::sort(ex.begin(), ex.end());
    std::sort(ey.begin(), ey.end());
    CHECK(ex == ey);
  }
}

bool generic_name(const std::string& s) {
  if (s == "anon" || s == "OBJ" || s == "RHS" || s == "RNG" || s == "BND" || s == "SOS") return true;
  if (s.size() >= 2 && std::strchr("xcsi", s[0])) {
    return std::all_of(s.begin() + 1, s.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); });
  }
  return false;
}

// No name of the original (other than ones that look generic themselves) is a token of the file.
void check_no_original_names(const Model& orig, const std::string& text) {
  std::set<std::string> tokens;
  std::istringstream in(text);
  for (std::string t; in >> t;) tokens.insert(t);
  std::vector<std::string> names = orig.col_names;
  names.insert(names.end(), orig.row_names.begin(), orig.row_names.end());
  names.push_back(orig.name);
  for (const SosConstraint& s : orig.sos) names.push_back(s.name);
  for (const IndicatorConstraint& c : orig.indicators) names.push_back(c.name);
  for (const std::string& nm : names) {
    if (nm.empty() || generic_name(nm)) continue;
    if (!CHECK(tokens.count(nm) == 0)) std::fprintf(stderr, "  original name '%s' in the anonymized file\n", nm.c_str());
  }
}

SolveResult solve_quietly(const Model& m) {
  SolverParams p;
  p.verbosity = 0;
  p.threads = 1;
  p.deterministic = true;
  p.time_limit = 60;
  return solve(m, p);
}

// Anonymize (shuffled), write, read back; check the copy exactly; solve both; map the copy's
// solution back to the original's columns.
void round_trip(const Model& orig, uint64_t seed, bool solve_both = true) {
  AnonMap map;
  AnonOptions opt;
  opt.seed = seed;
  const Model a = anonymize(orig, opt, &map);
  check_mapped(orig, a, map, false);
  Index inexact = -1;
  const std::string text = write_mps_string(a, &inexact);
  CHECK(inexact == 0);
  check_no_original_names(orig, text);
  std::vector<std::string> warnings;
  const Model b = read_mps_string(text, &warnings);
  CHECK(warnings.empty());
  check_mapped(orig, b, map, true);
  CHECK(col_multiset(b) == col_multiset(orig));
  CHECK(row_multiset(b) == row_multiset(orig));
  CHECK(b.name == "anon");
  for (Index j = 0; j < b.num_cols(); ++j) CHECK(b.col_names[j] == "x" + std::to_string(j + 1));
  if (!solve_both) return;
  const SolveResult ro = solve_quietly(orig), ra = solve_quietly(b);
  CHECK(ro.status == ra.status);
  if (ro.status == Status::Optimal && ra.status == Status::Optimal) {
    CHECK_NEAR(ra.objective, ro.objective, 1e-6 * std::max(1.0, std::abs(ro.objective)));
    // The copy's solution in the original's columns is a solution of the original with that objective.
    std::vector<double> x(static_cast<size_t>(orig.num_cols()));
    for (Index j = 0; j < b.num_cols(); ++j) x[static_cast<size_t>(map.col[j])] = ra.x[j];
    CHECK(orig.violation(x).max() <= 1e-5);
    CHECK_NEAR(orig.objective_value(x), ra.objective, 1e-9 * std::max(1.0, std::abs(ra.objective)));
  }
}

// Uniform in [-1, 1) with 53 random bits, the same on every standard library.
double signed_unit(std::mt19937_64& rng) { return static_cast<double>(rng() >> 10) * 0x1.0p-53 - 1.0; }

Model one_row(double lo, double hi) {
  ModelBuilder mb;
  mb.add_col("x", 1, 0, 1);
  const Index r = mb.add_row("r", lo, hi);
  mb.add_coef(r, 0, 1.0);
  return mb.build();
}

// hi, or the next double above it for which an MPS range gives [lo, hi] exactly (most pairs need
// no change; see write_mps on the ranged rows MPS cannot express).
double representable_upper(double lo, double hi) {
  if (!(lo > -kInf && hi < kInf && lo < hi)) return hi;
  for (int k = 0; k < 64; ++k, hi = std::nextafter(hi, kInf)) {
    Index inexact = -1;
    write_mps_string(one_row(lo, hi), &inexact);
    if (inexact == 0) return hi;
  }
  return hi;
}

// A feasible, bounded random MILP with full-precision data: rows built around a random point
// (inside the bounds, integral where needed) with random senses, ranges included. Names start with
// "secret" (none may reach the anonymized file).
Model random_milp(std::mt19937_64& rng) {
  std::uniform_real_distribution<double> U(0.0, 1.0);
  auto uni = [&](double a, double b) { return a + (b - a) * U(rng); };
  const Index n = 4 + static_cast<Index>(rng() % 10), m = 3 + static_cast<Index>(rng() % 8);
  ModelBuilder mb;
  std::vector<double> x0;
  for (Index j = 0; j < n; ++j) {
    const int t = static_cast<int>(rng() % 3);
    const std::string name = "secret_col_" + std::to_string(j);
    if (t == 0) {
      const double lb = uni(-5, 0), ub = lb + uni(0.5, 10);
      mb.add_col(name, uni(-10, 10), lb, ub, VarType::Continuous);
      x0.push_back(uni(lb, ub));
    } else if (t == 1) {
      const double lb = -static_cast<double>(rng() % 4), ub = lb + 1 + static_cast<double>(rng() % 5);
      mb.add_col(name, uni(-10, 10), lb, ub, VarType::Integer);
      x0.push_back(lb + static_cast<double>(rng() % static_cast<unsigned>(ub - lb + 1)));
    } else {
      mb.add_col(name, uni(-10, 10), 0, 1, VarType::Binary);
      x0.push_back(static_cast<double>(rng() % 2));
    }
  }
  for (Index i = 0; i < m; ++i) {
    std::vector<std::pair<Index, double>> e;
    double act = 0.0;
    for (Index j = 0; j < n; ++j) {
      if (U(rng) < 0.5) continue;
      const double v = uni(-5, 5);
      e.emplace_back(j, v);
      act += v * x0[j];
    }
    double lb = -kInf, ub = kInf;
    switch (rng() % 4) {
      case 0: ub = act + uni(0, 3); break;
      case 1: lb = act - uni(0, 3); break;
      case 2: lb = act - uni(0, 2), ub = representable_upper(lb, act + uni(0, 2)); break;
      default: lb = ub = act; break;
    }
    const Index r = mb.add_row("secret_row_" + std::to_string(i), lb, ub);
    for (const auto& [j, v] : e) mb.add_coef(r, j, v);
  }
  Model model = mb.build();
  model.name = "secret_model";
  model.obj_offset = uni(-100, 100);
  if (rng() % 2) model.sense = ObjSense::Maximize;
  return model;
}

Model read_instance(const std::string& rel) {
  std::vector<std::string> warnings;  // the lpsanity files test the reader's warnings
  return read_mps(std::string(MILP_SOURCE_DIR) + "/" + rel, &warnings);
}

}  // namespace

TEST(anon_random_milps_round_trip) {
  std::mt19937_64 rng(20261105);
  for (int t = 0; t < 25; ++t) round_trip(random_milp(rng), 1 + static_cast<uint64_t>(t));
}

TEST(anon_shuffles_and_is_deterministic) {
  std::mt19937_64 rng(7);
  Model m = random_milp(rng);
  while (m.num_cols() < 8 || m.num_rows() < 6) m = random_milp(rng);
  AnonOptions o1, o2;
  o2.seed = 2;
  AnonMap m1, m1b, m2;
  const std::string t1 = write_mps_string(anonymize(m, o1, &m1));
  CHECK(t1 == write_mps_string(anonymize(m, o1, &m1b)));  // same seed, same file
  CHECK(m1.col == m1b.col && m1.row == m1b.row);
  const std::string t2 = write_mps_string(anonymize(m, o2, &m2));
  CHECK(t1 != t2);
  CHECK(m1.col != m2.col || m1.row != m2.row);
  // Some column or row actually moved.
  bool moved = false;
  for (size_t k = 0; k < m1.col.size(); ++k) moved = moved || m1.col[k] != static_cast<Index>(k);
  for (size_t k = 0; k < m1.row.size(); ++k) moved = moved || m1.row[k] != static_cast<Index>(k);
  CHECK(moved);
}

TEST(anon_keep_order) {
  std::mt19937_64 rng(11);
  const Model m = random_milp(rng);
  AnonOptions o;
  o.keep_order = true;
  AnonMap map;
  const Model a = anonymize(m, o, &map);
  for (size_t k = 0; k < map.col.size(); ++k) CHECK(map.col[k] == static_cast<Index>(k));
  for (size_t k = 0; k < map.row.size(); ++k) CHECK(map.row[k] == static_cast<Index>(k));
  CHECK(a.obj == m.obj);
  CHECK(a.col_lb == m.col_lb);
  CHECK(a.col_ub == m.col_ub);
  CHECK(a.row_lb == m.row_lb);
  CHECK(a.row_ub == m.row_ub);
  CHECK(a.A.start == m.A.start);
  CHECK(a.A.index == m.A.index);
  CHECK(a.A.value == m.A.value);
  CHECK(a.row_names[0] == "c1");
  check_mapped(m, read_mps_string(write_mps_string(a)), map, true);
}

TEST(anon_map_inverts_permutation) {
  std::mt19937_64 rng(3);
  const Model m = random_milp(rng);
  AnonMap map;
  const Model a = anonymize(m, {}, &map);
  const std::string text = anon_map_string(m, map);
  std::istringstream in(text);
  std::string line;
  REQUIRE(std::getline(in, line));
  CHECK(line[0] == '#');
  // Columns then rows: "<new> <original>", the original name of the entry the map points to.
  for (Index j = 0; j < a.num_cols(); ++j) {
    REQUIRE(std::getline(in, line));
    CHECK(line == a.col_names[j] + " " + m.col_names[map.col[j]]);
  }
  for (Index i = 0; i < a.num_rows(); ++i) {
    REQUIRE(std::getline(in, line));
    CHECK(line == a.row_names[i] + " " + m.row_names[map.row[i]]);
  }
  CHECK(!std::getline(in, line));
  // Unnamed originals get their index.
  Model u = m;
  u.col_names[map.col[0]].clear();
  const std::string tu = anon_map_string(u, map);
  CHECK(tu.find("\nx1 C" + std::to_string(map.col[0]) + "\n") != std::string::npos);
}

// Exact numbers: values with full mantissas, huge and tiny magnitudes, infinite bounds and ranged
// rows survive the writer and the reader bit for bit. (Bounds and sides stay below 1e30, which MPS
// readers take as infinite.)
TEST(anon_numbers_exact) {
  std::mt19937_64 rng(5);
  auto wild = [&]() {
    const int e = static_cast<int>(rng() % 60) - 30;
    return std::ldexp(signed_unit(rng), e * 10);  // 2^-300 .. 2^290, full mantissas
  };
  auto wild_side = [&]() {
    const int e = static_cast<int>(rng() % 40) - 30;
    return std::ldexp(signed_unit(rng), e * 10);  // 2^-300 .. 2^90
  };
  ModelBuilder mb;
  const Index n = 30;
  for (Index j = 0; j < n; ++j) {
    double lb = wild_side(), ub = wild_side();
    if (lb > ub) std::swap(lb, ub);
    if (j % 7 == 0) lb = -kInf;
    if (j % 5 == 0) ub = kInf;
    mb.add_col("col" + std::to_string(j), wild(), lb, ub, VarType::Continuous);
  }
  mb.add_col("tiny", 4.9406564584124654e-324, 0.1, 0.30000000000000004);
  for (Index i = 0; i < 40; ++i) {
    double lb = wild_side(), ub = wild_side();
    if (lb > ub) std::swap(lb, ub);
    if (i % 4 == 1) lb = -kInf;
    if (i % 4 == 2) ub = kInf;
    if (i % 4 == 3) ub = lb;
    if (i % 4 == 0) ub = representable_upper(lb, ub);
    const Index r = mb.add_row("r" + std::to_string(i), lb, ub);
    for (Index j = 0; j < n; ++j)
      if (rng() % 3 == 0) mb.add_coef(r, j, wild());
  }
  const Index r = mb.add_row("decimal", -5.0, -1.8);  // -5 + (-1.8 + 5) is not -1.8: written as an L row
  mb.add_coef(r, 0, 1.0);
  Model m = mb.build();
  m.obj_offset = wild_side();
  round_trip(m, 9, false);
}

TEST(anon_lp_files) {
  for (const char* f : {"lpsanity_ranges_0", "lpsanity_objoffset_0", "lpsanity_objsense_maximize_0", "lpsanity_empty_col_0",
                        "lpsanity_extra_free_row_0", "lpsanity_freevars_0", "lpsanity_negative_ub_0",
                        "lpsanity_mip_mixed_bounds_0", "lpsanity_mip_marker_defaults_0", "lpsanity_infeasible_0"}) {
    round_trip(read_instance(std::string("benchmarks/instances/") + f + ".mps"), 1);
  }
}

TEST(anon_milp_files) {
  for (const char* f : {"knapsack_tiny_1", "randmilp_tiny_1", "setcover_tiny_1", "binpack_tiny_1", "gap_tiny_1"})
    round_trip(read_instance(std::string("benchmarks/instances/") + f + ".mps"), 3);
}

TEST(anon_quadratic) {
  for (const char* f : {"tests/data/maros_meszaros/HS21.mps", "tests/data/maros_meszaros/HS35.mps",
                        "tests/data/maros_meszaros/GENHS28.mps", "benchmarks/instances/portfolio_tiny_1.mps",
                        "benchmarks/instances/miqprand_tiny_1.mps"}) {
    const Model m = read_instance(f);
    REQUIRE(m.has_quadratic());
    round_trip(m, 4);
  }
}

TEST(anon_semicontinuous) {
  for (const char* f : {"semilot_tiny_1", "semilot_small_1"}) {
    const Model m = read_instance(std::string("benchmarks/instances/") + f + ".mps");
    REQUIRE(m.num_semi() > 0);
    round_trip(m, 5);
  }
}

TEST(anon_indicators_and_sos) {
  const Model ind = read_instance("benchmarks/instances/indsched_tiny_1.mps");
  REQUIRE(!ind.indicators.empty());
  round_trip(ind, 6);
  const Model pwl = read_instance("benchmarks/instances/pwlnet_tiny_1.mps");
  REQUIRE(!pwl.sos.empty() && !pwl.indicators.empty() && pwl.num_semi() > 0);
  round_trip(pwl, 7);
  // SOS2 with weights out of order and an indicator on value 0, built by hand.
  ModelBuilder mb;
  for (int j = 0; j < 6; ++j) mb.add_col("secret" + std::to_string(j), j % 2 ? -1.0 : -0.5, 0, 4);
  const Index z = mb.add_col("secret_z", 1, 0, 1, VarType::Binary);
  const Index r = mb.add_row("secret_sum", -kInf, 6);
  for (int j = 0; j < 6; ++j) mb.add_coef(r, j, 1.0 + 0.1 * j);
  Model m = mb.build();
  m.sos.push_back({"secret_sos", 2, {3, 1, 4, 0}, {3.5, 1.25, 4.0, 0.5}});
  m.sos.push_back({"secret_sos1", 1, {5, 2}, {2.0, 1.0}});
  m.indicators.push_back({"secret_ind", z, 0, {0, 5}, {1.0, 2.5}, -kInf, 3.0});
  m.indicators.push_back({"secret_ind2", z, 1, {1, 2}, {1.0, -1.0}, 0.5, kInf});
  round_trip(m, 8);
}
