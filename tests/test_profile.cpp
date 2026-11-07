// Tests for the model / run profile (solver/profile.hpp, popso --profile): the row classes on
// hand-built rows, the JSON (parsed by a strict minimal parser) for LPs, MILPs, QPs and logical
// models, no names in it, and the incumbent trace (which must not change the search).

#include <cctype>
#include <cmath>
#include <cstdlib>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "milp/io/mps_reader.hpp"
#include "milp/solver/profile.hpp"
#include "milp/solver/solver.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

// Strict JSON (RFC 8259) parser, enough to validate the profile and read fields from it.
struct Json {
  enum Type { Null, Bool, Number, String, Array, Object } type = Null;
  bool boolean = false;
  double number = 0.0;
  std::string str;
  std::vector<Json> arr;
  std::vector<std::pair<std::string, Json>> obj;

  const Json& operator[](const std::string& key) const {
    if (type != Object) throw std::runtime_error("not an object: " + key);
    for (const auto& [k, v] : obj)
      if (k == key) return v;
    throw std::runtime_error("missing key " + key);
  }
  bool has(const std::string& key) const {
    for (const auto& kv : obj)
      if (kv.first == key) return true;
    return false;
  }
};

class JsonParser {
 public:
  explicit JsonParser(const std::string& s) : s_(s) {}
  Json parse() {
    Json v = value();
    ws();
    if (i_ != s_.size()) fail("trailing characters");
    return v;
  }

 private:
  [[noreturn]] void fail(const char* what) { throw std::runtime_error(std::string("JSON: ") + what + " at " + std::to_string(i_)); }
  void ws() {
    while (i_ < s_.size() && (s_[i_] == ' ' || s_[i_] == '\t' || s_[i_] == '\n' || s_[i_] == '\r')) ++i_;
  }
  char peek() { return i_ < s_.size() ? s_[i_] : '\0'; }
  void expect(char c) {
    if (peek() != c) fail("unexpected character");
    ++i_;
  }
  bool literal(const char* w) {
    const std::string lit(w);
    if (s_.compare(i_, lit.size(), lit) != 0) return false;
    i_ += lit.size();
    return true;
  }
  Json value() {
    ws();
    Json v;
    const char c = peek();
    if (c == '{') {
      v.type = Json::Object;
      ++i_;
      ws();
      if (peek() == '}') {
        ++i_;
        return v;
      }
      for (;;) {
        ws();
        if (peek() != '"') fail("expected a key");
        std::string k = string();
        ws();
        expect(':');
        v.obj.emplace_back(std::move(k), value());
        ws();
        if (peek() == ',') {
          ++i_;
          continue;
        }
        expect('}');
        return v;
      }
    }
    if (c == '[') {
      v.type = Json::Array;
      ++i_;
      ws();
      if (peek() == ']') {
        ++i_;
        return v;
      }
      for (;;) {
        v.arr.push_back(value());
        ws();
        if (peek() == ',') {
          ++i_;
          continue;
        }
        expect(']');
        return v;
      }
    }
    if (c == '"') {
      v.type = Json::String;
      v.str = string();
      return v;
    }
    if (literal("true")) {
      v.type = Json::Bool;
      v.boolean = true;
      return v;
    }
    if (literal("false")) {
      v.type = Json::Bool;
      return v;
    }
    if (literal("null")) return v;
    v.type = Json::Number;
    v.number = number();
    return v;
  }
  std::string string() {
    expect('"');
    std::string out;
    for (;;) {
      if (i_ >= s_.size()) fail("unterminated string");
      const char c = s_[i_++];
      if (c == '"') return out;
      if (static_cast<unsigned char>(c) < 0x20) fail("control character in a string");
      if (c != '\\') {
        out += c;
        continue;
      }
      const char e = peek();
      ++i_;
      if (e == 'u') {
        if (i_ + 4 > s_.size()) fail("bad escape");
        for (int k = 0; k < 4; ++k)
          if (!std::isxdigit(static_cast<unsigned char>(s_[i_ + k]))) fail("bad escape");
        out += static_cast<char>(std::strtol(s_.substr(i_, 4).c_str(), nullptr, 16));
        i_ += 4;
      } else if (std::string("\"\\/bfnrt").find(e) != std::string::npos) {
        out += e;
      } else {
        fail("bad escape");
      }
    }
  }
  double number() {
    const size_t b = i_;
    if (peek() == '-') ++i_;
    if (peek() == '0') {
      ++i_;
    } else if (peek() >= '1' && peek() <= '9') {
      while (std::isdigit(static_cast<unsigned char>(peek()))) ++i_;
    } else {
      fail("bad value");
    }
    if (peek() == '.') {
      ++i_;
      if (!std::isdigit(static_cast<unsigned char>(peek()))) fail("bad fraction");
      while (std::isdigit(static_cast<unsigned char>(peek()))) ++i_;
    }
    if (peek() == 'e' || peek() == 'E') {
      ++i_;
      if (peek() == '+' || peek() == '-') ++i_;
      if (!std::isdigit(static_cast<unsigned char>(peek()))) fail("bad exponent");
      while (std::isdigit(static_cast<unsigned char>(peek()))) ++i_;
    }
    return std::strtod(s_.substr(b, i_ - b).c_str(), nullptr);
  }
  const std::string& s_;
  size_t i_ = 0;
};

// Every string (keys and values) of a JSON tree.
void strings_of(const Json& v, std::set<std::string>& out) {
  if (v.type == Json::String) out.insert(v.str);
  for (const Json& a : v.arr) strings_of(a, out);
  for (const auto& [k, x] : v.obj) {
    out.insert(k);
    strings_of(x, out);
  }
}

Model read_instance(const std::string& rel) {
  std::vector<std::string> warnings;
  return read_mps(std::string(MILP_SOURCE_DIR) + "/" + rel, &warnings);
}

SolverParams quiet_params() {
  SolverParams p;
  p.verbosity = 0;
  p.threads = 1;
  p.deterministic = true;
  p.trace_incumbents = true;
  p.time_limit = 60;
  return p;
}

// Solves, profiles and checks the JSON: it parses, has the schema's sections, agrees with the
// model and the result, and contains none of the model's names.
Json profile_of(const Model& m, const SolverParams& params, SolveResult* result_out = nullptr) {
  const SolveResult r = solve(m, params);
  const std::string text = profile_json(m, params, r, "test-version");
  Json js;
  try {
    js = JsonParser(text).parse();
  } catch (const std::exception& e) {
    CHECK(false);
    std::fprintf(stderr, "  %s\n%s\n", e.what(), text.c_str());
    throw;
  }
  CHECK(js["profile_version"].number == 1);
  CHECK(js["solver"]["version"].str == "test-version");
  const Json& model = js["model"];
  CHECK(model["rows"].number == m.num_rows());
  CHECK(model["cols"].number == m.num_cols());
  CHECK(model["nonzeros"].number == m.num_nonzeros());
  double classes = 0.0, types = 0.0;
  for (const auto& kv : model["row_classes"].obj) classes += kv.second.number;
  for (const auto& kv : model["col_types"].obj) types += kv.second.number;
  CHECK(classes == m.num_rows());
  CHECK(types == m.num_cols());
  CHECK(model["row_classes"].obj.size() == static_cast<size_t>(kNumRowClasses));
  for (const char* sec : {"numerics", "presolve", "run"}) CHECK(js.has(sec));
  const Json& run = js["run"];
  CHECK(run["status"].str == to_string(r.status));
  CHECK(run["has_solution"].boolean == r.has_solution);
  if (r.has_solution) CHECK(run["objective"].number == r.objective);
  CHECK(run["nodes"].number == static_cast<double>(r.stats.nodes));
  CHECK(run["incumbents"].arr.size() == r.incumbent_trace.size());
  // No names: none of the model's names is a key or a string value, and none occurs in the text
  // (names that could occur by chance, like "x1", are not used by these models' checks below).
  std::set<std::string> strings;
  strings_of(js, strings);
  std::vector<std::string> names = m.col_names;
  names.insert(names.end(), m.row_names.begin(), m.row_names.end());
  names.push_back(m.name);
  for (const SosConstraint& s : m.sos) names.push_back(s.name);
  for (const IndicatorConstraint& c : m.indicators) names.push_back(c.name);
  for (const std::string& nm : names) {
    if (!nm.empty() && !CHECK(strings.count(nm) == 0)) std::fprintf(stderr, "  name '%s' in the profile\n", nm.c_str());
  }
  if (result_out) *result_out = r;
  return js;
}

}  // namespace

// One row of each class, in the table's order, plus the edge cases the rules spell out.
TEST(profile_row_classes) {
  ModelBuilder mb;
  const Index b1 = mb.add_col("b1", 0, 0, 1, VarType::Binary), b2 = mb.add_col("b2", 0, 0, 1, VarType::Binary);
  const Index b3 = mb.add_col("b3", 0, 0, 1, VarType::Binary), b4 = mb.add_col("b4", 0, 0, 1, VarType::Integer);
  const Index i1 = mb.add_col("i1", 0, 0, 10, VarType::Integer), i2 = mb.add_col("i2", 0, 0, 10, VarType::Integer);
  const Index c1 = mb.add_col("c1", 0, 0, 10), c2 = mb.add_col("c2", 0, 0, 10);
  const Index sc = mb.add_col("sc", 0, 1, 10, VarType::SemiContinuous);
  const Index si = mb.add_col("si", 0, 0, 1, VarType::SemiInteger);
  std::vector<RowClass> expect;
  auto row = [&](double lb, double ub, std::vector<std::pair<Index, double>> e, RowClass c) {
    const Index r = mb.add_row("r" + std::to_string(expect.size()), lb, ub);
    for (const auto& [j, v] : e) mb.add_coef(r, j, v);
    expect.push_back(c);
  };
  row(0, 1, {}, RowClass::Empty);
  row(-kInf, kInf, {{b1, 1}, {b2, 1}}, RowClass::Free);
  row(-kInf, 3, {{c1, 2}}, RowClass::Singleton);
  row(1, 1, {{c1, 1}, {c2, 2}}, RowClass::Aggregation);
  row(-kInf, 0, {{i1, 1}, {i2, -1}}, RowClass::Precedence);
  row(-kInf, 0, {{c1, 1}, {b1, -10}}, RowClass::VariableBound);
  row(1, 1, {{b1, 1}, {b2, 1}, {b3, 1}}, RowClass::SetPartitioning);
  row(-kInf, 1, {{b1, 1}, {b2, 1}, {b3, 1}}, RowClass::SetPacking);
  row(1, kInf, {{b1, 1}, {b2, 1}, {b3, 1}}, RowClass::SetCovering);
  row(2, 2, {{b1, 1}, {b2, 1}, {b3, 1}, {b4, 1}}, RowClass::Cardinality);
  row(-kInf, 2, {{b1, 1}, {b2, 1}, {b3, 1}, {b4, 1}}, RowClass::InvariantKnapsack);
  row(5, 5, {{b1, 2}, {b2, 3}, {b3, 4}}, RowClass::EquationKnapsack);
  row(-kInf, 5, {{b1, 2}, {b2, 3}, {b3, 5}}, RowClass::BinPacking);
  row(-kInf, 6, {{b1, 2}, {b2, 3}, {b3, 4}}, RowClass::Knapsack);
  row(-kInf, 7, {{i1, 2}, {i2, 3}, {b1, 1}}, RowClass::IntegerKnapsack);
  row(-kInf, 3, {{b1, 1.5}, {c1, 1}, {c2, -1}}, RowClass::MixedBinary);
  row(-kInf, 4, {{i1, 1}, {c1, 1}, {b1, 0.5}}, RowClass::GeneralLinear);
  // Edge cases.
  row(-kInf, -1, {{b1, 1}, {b2, -1}, {b3, -1}}, RowClass::SetPacking);  // b1 + (1-b2) + (1-b3) <= 1
  row(2, kInf, {{b1, 1}, {b2, 1}, {b3, 1}}, RowClass::SetPacking);      // at most one of them is 0
  row(2, kInf, {{b1, 1}, {b2, 1}, {b3, 1}, {b4, 1}, {i1, 0}}, RowClass::InvariantKnapsack);  // explicit zero ignored
  row(1, 2, {{b1, 1}, {b2, 1}, {b3, 1}}, RowClass::SetCovering);        // ranged: the earlier class of its sides
  row(-kInf, 1, {{b1, 1}, {b2, 1}}, RowClass::VariableBound);           // two entries: before set packing
  row(-kInf, 0, {{b1, 1}, {b2, -1}}, RowClass::Precedence);             // ... and precedence before that
  row(-kInf, 0, {{c1, 1}, {i1, -1}}, RowClass::GeneralLinear);       // kinds differ, no binary: general
  row(-kInf, 1, {{c1, 1}, {c2, 1}, {sc, 1}}, RowClass::MixedBinary);    // continuous only (semi-continuous too)
  row(-kInf, 4, {{si, 1}, {i1, 1}, {i2, 2}}, RowClass::IntegerKnapsack);  // semi-integer counts as integer
  row(-kInf, 4, {{si, 1}, {b1, 1}, {b2, 1}}, RowClass::IntegerKnapsack);  // ... even with bounds [0, 1]
  row(-kInf, 1.5, {{b1, 3}, {b2, 2}, {b3, 4}}, RowClass::MixedBinary);  // fractional rhs: no knapsack
  row(0, 0, {{b1, 1}, {b2, 1}, {b3, 1}}, RowClass::MixedBinary);        // b' = 0: no cardinality (k >= 2)
  const Model m = mb.build();
  const std::vector<RowClass> got = classify_rows(m);
  REQUIRE(got.size() == expect.size());
  for (size_t i = 0; i < got.size(); ++i) {
    if (!CHECK(got[i] == expect[i]))
      std::fprintf(stderr, "  row %zu: %s, expected %s\n", i, to_string(got[i]), to_string(expect[i]));
  }
  // Each class appears in the JSON counts.
  const Json js = JsonParser(profile_json(m, quiet_params(), SolveResult{}, "v")).parse();
  for (int c = 0; c < kNumRowClasses; ++c) CHECK(js["model"]["row_classes"][to_string(static_cast<RowClass>(c))].number >= 1);
}

TEST(profile_json_lp_milp_qp_logical) {
  // LP with ranged rows.
  {
    const Json js = profile_of(read_instance("benchmarks/instances/lpsanity_ranges_0.mps"), quiet_params());
    CHECK(js["run"]["status"].str == "OPTIMAL");
    CHECK(js["model"]["row_senses"]["ranged"].number >= 1);
  }
  // MILP: presolve sizes, incumbents with sources.
  {
    SolveResult r;
    const Json js = profile_of(read_instance("benchmarks/instances/knapsack_tiny_1.mps"), quiet_params(), &r);
    CHECK(js["presolve"]["ran"].boolean);
    CHECK(js["presolve"]["rows"].number <= js["presolve"]["rows_before"].number);
    CHECK(js["presolve"]["rows_before"].number == r.stats.presolve.rows_before);
    REQUIRE(!js["run"]["incumbents"].arr.empty());
    for (const Json& e : js["run"]["incumbents"].arr) {
      CHECK(!e["source"].str.empty());
      CHECK(e["count"].number >= 1);
      CHECK(e["last_time"].number >= e["time"].number);
    }
    CHECK(js["model"]["row_classes"]["knapsack"].number + js["model"]["row_classes"]["bin_packing"].number +
              js["model"]["row_classes"]["invariant_knapsack"].number >= 1);
  }
  // Continuous QP and MIQP.
  {
    const Json js = profile_of(read_instance("tests/data/maros_meszaros/HS21.mps"), quiet_params());
    CHECK(js["model"]["quadratic"]["nonzeros"].number >= 1);
    CHECK(js["run"]["status"].str == "OPTIMAL");
    CHECK(js["run"]["incumbents"].arr.size() == 1);  // the solution returned ("final")
  }
  {
    const Json js = profile_of(read_instance("benchmarks/instances/miqprand_tiny_1.mps"), quiet_params());
    CHECK(js["model"]["quadratic"]["diagonal"].number >= 1);
    CHECK(!js["run"]["incumbents"].arr.empty());
  }
  // SOS, indicators, semi-continuous columns.
  {
    const Model m = read_instance("benchmarks/instances/pwlnet_tiny_1.mps");
    const Json js = profile_of(m, quiet_params());
    CHECK(js["model"]["sos1"].number + js["model"]["sos2"].number == static_cast<double>(m.sos.size()));
    CHECK(js["model"]["indicators"].number == static_cast<double>(m.indicators.size()));
    CHECK(js["model"]["col_types"]["semi_continuous"].number == m.num_semi());
  }
}

// Distinctive names (and characters JSON must escape in them) never reach the profile.
TEST(profile_has_no_names) {
  ModelBuilder mb;
  const Index x = mb.add_col("secret\"col\\1", -1, 0, 4, VarType::Integer);
  const Index y = mb.add_col("secret_col_2", -2, 0, 3, VarType::Integer);
  const Index r = mb.add_row("secret_row\n1", -kInf, 5.5);
  mb.add_coef(r, x, 1.0);
  mb.add_coef(r, y, 1.5);
  Model m = mb.build();
  m.name = "secret model";
  SolveResult res = solve(m, quiet_params());
  const std::string text = profile_json(m, quiet_params(), res, "v\"1");
  CHECK(text.find("secret") == std::string::npos);
  const Json js = JsonParser(text).parse();  // the version's quote is escaped
  CHECK(js["solver"]["version"].str == "v\"1");
  profile_of(m, quiet_params());
}

// The trace: improving objectives in the model's sense, nondecreasing times, ending at the result.
// Tracing only records: the deterministic search is the same with and without it.
TEST(profile_incumbent_trace) {
  for (const char* f : {"benchmarks/instances/gap_tiny_1.mps", "benchmarks/instances/setcover_tiny_1.mps"}) {
    Model m = read_instance(f);
    for (int sense = 0; sense < 2; ++sense) {
      if (sense == 1) {
        // The same problem as a maximization of the negated objective.
        m.sense = ObjSense::Maximize;
        for (double& c : m.obj) c = -c;
      }
      SolverParams p = quiet_params();
      p.threads = 2;
      const SolveResult r = solve(m, p);
      REQUIRE(r.status == Status::Optimal);
      REQUIRE(!r.incumbent_trace.empty());
      const double s = m.sense == ObjSense::Maximize ? -1.0 : 1.0;
      double prev_t = 0.0, prev_obj = kInf;
      for (const IncumbentEvent& e : r.incumbent_trace) {
        CHECK(e.time >= prev_t);
        CHECK(e.last_time >= e.time);
        CHECK(s * e.objective < prev_obj);
        CHECK(s * e.last_objective <= s * e.objective);
        CHECK(e.count >= 1);
        prev_t = e.last_time;
        prev_obj = s * e.last_objective;
      }
      CHECK_NEAR(r.incumbent_trace.back().last_objective, r.objective, 1e-6 * std::max(1.0, std::abs(r.objective)));
      SolverParams q = p;
      q.trace_incumbents = false;
      const SolveResult r2 = solve(m, q);
      CHECK(r2.incumbent_trace.empty());
      CHECK(r2.stats.nodes == r.stats.nodes);
      CHECK(r2.stats.lp_iterations == r.stats.lp_iterations);
      CHECK(r2.objective == r.objective);
    }
  }
}

// A solve stopped by the time limit still gives a complete profile.
TEST(profile_time_limit) {
  const Model m = read_instance("benchmarks/instances/cflp_medium_1.mps");
  SolverParams p = quiet_params();
  p.time_limit = 0.05;
  const Json js = profile_of(m, p);
  CHECK(js["run"]["status"].str == "TIME_LIMIT");
  CHECK(js["solver"]["time_limit"].number == 0.05);
}
