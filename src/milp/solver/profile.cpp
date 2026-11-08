#include "milp/solver/profile.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <thread>

#include "milp/io/charconv_compat.hpp"

namespace milp {

const char* to_string(RowClass c) {
  switch (c) {
    case RowClass::Empty: return "empty";
    case RowClass::Free: return "free";
    case RowClass::Singleton: return "singleton";
    case RowClass::Aggregation: return "aggregation";
    case RowClass::Precedence: return "precedence";
    case RowClass::VariableBound: return "variable_bound";
    case RowClass::SetPartitioning: return "set_partitioning";
    case RowClass::SetPacking: return "set_packing";
    case RowClass::SetCovering: return "set_covering";
    case RowClass::Cardinality: return "cardinality";
    case RowClass::InvariantKnapsack: return "invariant_knapsack";
    case RowClass::EquationKnapsack: return "equation_knapsack";
    case RowClass::BinPacking: return "bin_packing";
    case RowClass::Knapsack: return "knapsack";
    case RowClass::IntegerKnapsack: return "integer_knapsack";
    case RowClass::MixedBinary: return "mixed_binary";
    case RowClass::GeneralLinear: return "general_linear";
  }
  return "general_linear";
}

namespace {

enum class Kind : uint8_t { Binary, Integer, Continuous };

Kind column_kind(const Model& m, Index j) {
  const VarType t = m.col_type[j];
  if (t == VarType::Continuous || t == VarType::SemiContinuous) return Kind::Continuous;
  if (t != VarType::SemiInteger && m.col_lb[j] > -1.0 && m.col_ub[j] < 2.0) return Kind::Binary;
  return Kind::Integer;
}

bool near(double a, double b) { return std::abs(a - b) <= 1e-9 * std::max({1.0, std::abs(a), std::abs(b)}); }
bool integral(double v) { return std::isfinite(v) && near(v, std::round(v)); }
bool integral_at_least(double v, double k) { return integral(v) && std::round(v) >= k; }

// One row: its entries (nonzeros) and what the rules need about them.
struct RowView {
  std::vector<double> a;
  std::vector<Kind> kind;
  std::vector<double> col_lb;
  double lb = -kInf, ub = kInf;
};

// Classes 8, 11, 13-15 for the <= form sign * a^T x <= beta; GeneralLinear if none matches.
RowClass le_form_class(const RowView& r, double sign, double beta, bool all_bin, bool same_abs, bool int_coefs) {
  const size_t n = r.a.size();
  if (all_bin) {
    double neg_sum = 0.0;  // complementing the binaries with negative coefficients
    for (double v : r.a)
      if (sign * v < 0.0) neg_sum += std::abs(v);
    const double b = beta + neg_sum;
    if (same_abs) {
      const double s = std::abs(r.a[0]);
      if (near(b / s, 1.0)) return RowClass::SetPacking;
      if (integral_at_least(b / s, 2.0)) return RowClass::InvariantKnapsack;
    }
    if (int_coefs && integral_at_least(b, 2.0)) {
      for (double v : r.a)
        if (near(std::abs(v), b)) return RowClass::BinPacking;
      return RowClass::Knapsack;
    }
  }
  bool int_knapsack = integral_at_least(beta, 0.0);
  for (size_t k = 0; k < n && int_knapsack; ++k) {
    const double c = sign * r.a[k];
    int_knapsack = r.kind[k] != Kind::Continuous && r.col_lb[k] >= 0.0 && c > 0.0 && integral(c);
  }
  if (int_knapsack) return RowClass::IntegerKnapsack;
  return RowClass::GeneralLinear;
}

RowClass classify(const RowView& r) {
  const size_t n = r.a.size();
  if (n == 0) return RowClass::Empty;
  const bool has_lb = r.lb > -kInf, has_ub = r.ub < kInf;
  if (!has_lb && !has_ub) return RowClass::Free;
  if (n == 1) return RowClass::Singleton;
  const bool equation = r.lb == r.ub;
  if (n == 2 && equation) return RowClass::Aggregation;
  if (n == 2) {
    if (r.kind[0] == r.kind[1] && near(r.a[0], -r.a[1])) return RowClass::Precedence;
    if (r.kind[0] == Kind::Binary || r.kind[1] == Kind::Binary) return RowClass::VariableBound;
  }
  bool all_bin = true, same_abs = true, int_coefs = true, bin_or_cont = true;
  int nneg = 0;
  const double s = std::abs(r.a[0]);
  for (size_t k = 0; k < n; ++k) {
    all_bin = all_bin && r.kind[k] == Kind::Binary;
    bin_or_cont = bin_or_cont && r.kind[k] != Kind::Integer;
    same_abs = same_abs && near(std::abs(r.a[k]), s);
    int_coefs = int_coefs && integral(r.a[k]);
    nneg += r.a[k] < 0.0;
  }
  const RowClass rest = bin_or_cont ? RowClass::MixedBinary : RowClass::GeneralLinear;
  if (equation) {
    if (all_bin) {
      double neg_sum = 0.0;
      for (double v : r.a)
        if (v < 0.0) neg_sum += std::abs(v);
      const double b = r.ub + neg_sum;
      if (same_abs && near(b / s, 1.0)) return RowClass::SetPartitioning;
      if (same_abs && integral_at_least(b / s, 2.0)) return RowClass::Cardinality;
      if (int_coefs && integral_at_least(b, 2.0)) return RowClass::EquationKnapsack;
    }
    return rest;
  }
  RowClass best = rest;
  if (has_ub) best = std::min(best, le_form_class(r, 1.0, r.ub, all_bin, same_abs, int_coefs));
  if (has_lb) {
    best = std::min(best, le_form_class(r, -1.0, -r.lb, all_bin, same_abs, int_coefs));
    if (all_bin && same_abs && near((r.lb + s * nneg) / s, 1.0)) best = std::min(best, RowClass::SetCovering);
  }
  return best;
}

// Minimal JSON writer: objects and arrays, two-space indentation, compact objects inside arrays.
class JsonWriter {
 public:
  JsonWriter() {
    out_ += '{';
    stack_.push_back({true, false});
  }
  void begin(const char* key, char open, bool compact = false) {
    item(key);
    out_ += open;
    stack_.push_back({true, compact});
  }
  void end(char close) {
    const Level lv = stack_.back();
    stack_.pop_back();
    if (!lv.first && !lv.compact) newline();
    out_ += close;
  }
  void num(const char* key, double v) {
    item(key);
    if (!std::isfinite(v)) {
      out_ += "null";
      return;
    }
    char buf[40];
    out_.append(buf, static_cast<size_t>(io_detail::format_double(buf, v)));
  }
  void integer(const char* key, long long v) {
    item(key);
    out_ += std::to_string(v);
  }
  void boolean(const char* key, bool v) {
    item(key);
    out_ += v ? "true" : "false";
  }
  void str(const char* key, const std::string& v) {
    item(key);
    quote(v);
  }
  std::string finish() {
    end('}');
    out_ += '\n';
    return std::move(out_);
  }

 private:
  struct Level {
    bool first, compact;
  };
  void newline() {
    out_ += '\n';
    out_.append(2 * stack_.size(), ' ');
  }
  void item(const char* key) {
    Level& lv = stack_.back();
    if (!lv.first) out_ += lv.compact ? ", " : ",";
    if (!lv.compact) newline();
    lv.first = false;
    if (key) {
      quote(key);
      out_ += ": ";
    }
  }
  void quote(const std::string& s) {
    out_ += '"';
    for (const char ch : s) {
      const unsigned char c = static_cast<unsigned char>(ch);
      if (c == '"' || c == '\\') {
        out_ += '\\';
        out_ += ch;
      } else if (c < 0x20) {
        char buf[8];
        std::snprintf(buf, sizeof(buf), "\\u%04x", c);
        out_ += buf;
      } else {
        out_ += ch;
      }
    }
    out_ += '"';
  }
  std::string out_;
  std::vector<Level> stack_;
};

// Smallest and largest of the nonzero absolute values fed to it.
struct AbsRange {
  double lo = kInf, hi = 0.0;
  void add(double v) {
    v = std::abs(v);
    if (v == 0.0 || !std::isfinite(v)) return;
    lo = std::min(lo, v);
    hi = std::max(hi, v);
  }
  void write(JsonWriter& js, const char* key) const {
    js.begin(key, '{', true);
    js.num("min", hi > 0.0 ? lo : kInf);
    js.num("max", hi > 0.0 ? hi : kInf);
    js.end('}');
  }
};

void write_lengths(JsonWriter& js, const char* key, std::vector<Index> len) {
  js.begin(key, '{', true);
  if (len.empty()) {
    js.num("min", kInf), js.num("median", kInf), js.num("max", kInf);
  } else {
    std::sort(len.begin(), len.end());
    const size_t k = len.size();
    js.integer("min", len.front());
    js.num("median", k % 2 ? static_cast<double>(len[k / 2]) : 0.5 * (static_cast<double>(len[k / 2 - 1]) + len[k / 2]));
    js.integer("max", len.back());
  }
  js.end('}');
}

}  // namespace

std::vector<RowClass> classify_rows(const Model& model) {
  const SparseMatrix R = model.A.transpose();
  std::vector<Kind> kind(static_cast<size_t>(model.num_cols()));
  for (Index j = 0; j < model.num_cols(); ++j) kind[j] = column_kind(model, j);
  std::vector<RowClass> out(static_cast<size_t>(model.num_rows()));
  RowView r;
  for (Index i = 0; i < model.num_rows(); ++i) {
    r.a.clear(), r.kind.clear(), r.col_lb.clear();
    for (Index p = R.begin(i); p < R.end(i); ++p) {
      if (R.value[p] == 0.0) continue;
      const Index j = R.index[p];
      r.a.push_back(R.value[p]);
      r.kind.push_back(kind[j]);
      r.col_lb.push_back(model.col_lb[j]);
    }
    r.lb = model.row_lb[i];
    r.ub = model.row_ub[i];
    out[i] = classify(r);
  }
  return out;
}

std::string profile_json(const Model& model, const SolverParams& params, const SolveResult& result,
                         const std::string& version) {
  const Index n = model.num_cols(), m = model.num_rows();
  JsonWriter js;
  js.integer("profile_version", 1);

  js.begin("solver", '{');
  js.str("name", "popso");
  js.str("version", version);
  const long long threads = params.threads > 0 ? params.threads : std::thread::hardware_concurrency();
  js.integer("threads", std::max(1LL, threads));
  js.num("time_limit", params.time_limit);
  js.boolean("deterministic", params.deterministic);
  js.integer("seed", params.seed);
  js.integer("racing", params.racing);
  js.boolean("presolve", params.presolve);
  js.num("mip_rel_gap", params.mip_rel_gap);
  js.str("lp_method", to_string(params.lp_method));
  js.end('}');

  // Columns.
  long long types[5] = {0, 0, 0, 0, 0};  // continuous, binary, integer, semi-continuous, semi-integer
  long long bounds[5] = {0, 0, 0, 0, 0};  // free, lower, upper, boxed, fixed
  AbsRange obj_range, bound_range, matrix_range, rhs_range, q_range;
  long long obj_nz = 0;
  std::vector<Index> col_len(static_cast<size_t>(n));
  for (Index j = 0; j < n; ++j) {
    const VarType t = model.col_type[j];
    if (t == VarType::SemiContinuous) ++types[3];
    else if (t == VarType::SemiInteger) ++types[4];
    else if (t == VarType::Continuous) ++types[0];
    else ++types[column_kind(model, j) == Kind::Binary ? 1 : 2];
    const double lo = model.col_lb[j], hi = model.col_ub[j];
    if (lo == hi) ++bounds[4];
    else if (lo > -kInf && hi < kInf) ++bounds[3];
    else if (lo > -kInf) ++bounds[1];
    else if (hi < kInf) ++bounds[2];
    else ++bounds[0];
    bound_range.add(lo);
    bound_range.add(hi);
    obj_range.add(model.obj[j]);
    obj_nz += model.obj[j] != 0.0;
    col_len[j] = model.A.end(j) - model.A.begin(j);
    for (Index p = model.A.begin(j); p < model.A.end(j); ++p) matrix_range.add(model.A.value[p]);
  }

  // Rows.
  long long senses[5] = {0, 0, 0, 0, 0};  // equal, less, greater, ranged, free
  for (Index i = 0; i < m; ++i) {
    const double lo = model.row_lb[i], hi = model.row_ub[i];
    if (lo == hi) ++senses[0];
    else if (lo > -kInf && hi < kInf) ++senses[3];
    else if (hi < kInf) ++senses[1];
    else if (lo > -kInf) ++senses[2];
    else ++senses[4];
    rhs_range.add(lo);
    rhs_range.add(hi);
  }
  long long classes[kNumRowClasses] = {};
  for (RowClass c : classify_rows(model)) ++classes[static_cast<int>(c)];
  const SparseMatrix R = model.A.transpose();
  std::vector<Index> row_len(static_cast<size_t>(m));
  long long wide_range = 0, tiny_entries = 0, continuous_only = 0;
  for (Index i = 0; i < m; ++i) {
    row_len[i] = R.end(i) - R.begin(i);
    double lo = kInf, hi = 0.0;
    bool cont = R.begin(i) < R.end(i);
    for (Index p = R.begin(i); p < R.end(i); ++p) {
      const double v = std::abs(R.value[p]);
      if (v > 0.0) lo = std::min(lo, v), hi = std::max(hi, v);
      cont = cont && column_kind(model, R.index[p]) == Kind::Continuous;
    }
    if (hi > 0.0 && hi > 1e6 * lo) ++wide_range;
    if (hi > 0.0 && lo < 1e-9 * hi) ++tiny_entries;
    continuous_only += cont;
  }

  // Logical constraints and the quadratic objective.
  long long sos_count[2] = {0, 0}, sos_members = 0, ind_nz = 0;
  for (const SosConstraint& s : model.sos) {
    ++sos_count[s.type == 2 ? 1 : 0];
    sos_members += static_cast<long long>(s.cols.size());
  }
  for (const IndicatorConstraint& c : model.indicators) ind_nz += static_cast<long long>(c.cols.size());
  const SparseMatrix Q = model.q_lower();
  long long q_diag = 0, q_cols = 0;
  {
    std::vector<char> in_q(static_cast<size_t>(n), 0);
    for (Index j = 0; j < Q.num_major; ++j) {
      for (Index p = Q.begin(j); p < Q.end(j); ++p) {
        q_diag += Q.index[p] == j;
        in_q[j] = in_q[Q.index[p]] = 1;
        q_range.add(Q.value[p]);
      }
    }
    for (char c : in_q) q_cols += c;
  }
  // Quadratic constraint terms (Model::qr_*), merged per (row, i, j).
  long long qr_rows = 0, qr_terms = 0, qr_diag = 0, qr_cols = 0, qr_max = 0;
  AbsRange qr_range;
  {
    Model qm;
    qm.qr_row = model.qr_row, qm.qr_i = model.qr_i, qm.qr_j = model.qr_j, qm.qr_val = model.qr_val;
    qm.compact_qrows();
    std::vector<char> in_q(static_cast<size_t>(n), 0);
    long long run = 0;
    for (size_t k = 0; k < qm.qr_val.size(); ++k) {
      if (k == 0 || qm.qr_row[k] != qm.qr_row[k - 1]) ++qr_rows, run = 0;
      qr_max = std::max(qr_max, ++run);
      ++qr_terms;
      qr_diag += qm.qr_i[k] == qm.qr_j[k];
      in_q[qm.qr_i[k]] = in_q[qm.qr_j[k]] = 1;
      qr_range.add(qm.qr_val[k]);
    }
    for (char c : in_q) qr_cols += c;
  }

  js.begin("model", '{');
  js.str("sense", model.sense == ObjSense::Maximize ? "maximize" : "minimize");
  js.integer("rows", m);
  js.integer("cols", n);
  js.integer("nonzeros", model.num_nonzeros());
  js.integer("objective_nonzeros", obj_nz);
  js.begin("col_types", '{', true);
  js.integer("continuous", types[0]);
  js.integer("binary", types[1]);
  js.integer("integer", types[2]);
  js.integer("semi_continuous", types[3]);
  js.integer("semi_integer", types[4]);
  js.end('}');
  js.begin("col_bounds", '{', true);
  js.integer("free", bounds[0]);
  js.integer("lower", bounds[1]);
  js.integer("upper", bounds[2]);
  js.integer("boxed", bounds[3]);
  js.integer("fixed", bounds[4]);
  js.end('}');
  js.begin("row_senses", '{', true);
  js.integer("equal", senses[0]);
  js.integer("less", senses[1]);
  js.integer("greater", senses[2]);
  js.integer("ranged", senses[3]);
  js.integer("free", senses[4]);
  js.end('}');
  js.begin("row_classes", '{');
  for (int c = 0; c < kNumRowClasses; ++c) js.integer(to_string(static_cast<RowClass>(c)), classes[c]);
  js.end('}');
  js.integer("rows_continuous_only", continuous_only);
  js.integer("sos1", sos_count[0]);
  js.integer("sos2", sos_count[1]);
  js.integer("sos_members", sos_members);
  js.integer("indicators", static_cast<long long>(model.indicators.size()));
  js.integer("indicator_nonzeros", ind_nz);
  js.begin("quadratic", '{', true);
  js.integer("nonzeros", Q.nnz());
  js.integer("diagonal", q_diag);
  js.integer("cols", q_cols);
  js.end('}');
  js.begin("quadratic_rows", '{', true);
  js.integer("rows", qr_rows);
  js.integer("terms", qr_terms);
  js.integer("diagonal", qr_diag);
  js.integer("cols", qr_cols);
  js.integer("max_row_terms", qr_max);
  js.end('}');
  js.end('}');

  js.begin("numerics", '{');
  matrix_range.write(js, "matrix");
  obj_range.write(js, "objective");
  bound_range.write(js, "bounds");
  rhs_range.write(js, "rhs");
  q_range.write(js, "quadratic");
  qr_range.write(js, "quadratic_rows");
  js.integer("rows_coef_range_above_1e6", wide_range);
  js.integer("rows_with_tiny_entries", tiny_entries);
  write_lengths(js, "row_length", std::move(row_len));
  write_lengths(js, "col_length", std::move(col_len));
  js.end('}');

  const PresolveSummary& ps = result.stats.presolve;
  js.begin("presolve", '{');
  js.boolean("ran", ps.ran);
  if (ps.ran) {
    js.integer("rows_before", ps.rows_before);
    js.integer("cols_before", ps.cols_before);
    js.integer("nonzeros_before", ps.nnz_before);
    js.integer("rows", ps.rows);
    js.integer("cols", ps.cols);
    js.integer("nonzeros", ps.nnz);
    js.num("time", ps.time);
  }
  js.end('}');

  js.begin("run", '{');
  js.str("status", to_string(result.status));
  js.boolean("has_solution", result.has_solution);
  js.num("objective", result.has_solution ? result.objective : kInf);
  js.num("bound", result.bound);
  double gap = kInf;
  if (result.has_solution && std::isfinite(result.bound)) {
    const double diff = std::abs(result.objective - result.bound);
    gap = diff == 0.0 ? 0.0 : diff / std::max(std::abs(result.objective), 1e-10);
  }
  js.num("gap", gap);
  js.num("time", result.stats.time);
  js.integer("nodes", result.stats.nodes);
  js.integer("lp_iterations", result.stats.lp_iterations);
  js.integer("barrier_iterations", result.stats.barrier_iterations);
  js.num("root_bound", result.stats.root_bound);
  js.integer("restarts", result.stats.restarts);
  js.integer("tree_restarts", result.stats.tree_restarts);
  js.begin("incumbents", '[');
  for (const IncumbentEvent& e : result.incumbent_trace) {
    js.begin(nullptr, '{', true);
    js.num("time", e.time);
    js.num("objective", e.objective);
    js.str("source", e.source);
    js.integer("count", e.count);
    js.num("last_time", e.last_time);
    js.num("last_objective", e.last_objective);
    js.end('}');
  }
  js.end(']');
  js.end('}');
  return js.finish();
}

void write_profile(const std::string& path, const Model& model, const SolverParams& params, const SolveResult& result,
                   const std::string& version) {
  const std::string text = profile_json(model, params, result, version);
  std::FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) throw std::runtime_error(path + ": cannot open file for writing");
  const bool wrote = std::fwrite(text.data(), 1, text.size(), f) == text.size();
  const bool closed = std::fclose(f) == 0;
  if (!wrote || !closed) throw std::runtime_error(path + ": write error");
}

}  // namespace milp
