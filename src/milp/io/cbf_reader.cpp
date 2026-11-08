#include "milp/io/cbf_reader.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <string_view>
#include <tuple>
#include <utility>
#include <vector>

#include "milp/io/gzip.hpp"

namespace milp {

namespace {

struct ConeChunk {
  std::string kind;
  Index size = 0;
};

class CbfParser {
 public:
  CbfParser(const std::string& text, std::string source) : source_(std::move(source)) {
    size_t line_no = 0;
    for (size_t p = 0; p <= text.size();) {
      size_t e = text.find('\n', p);
      if (e == std::string::npos) e = text.size();
      ++line_no;
      std::string_view l(text.data() + p, e - p);
      while (!l.empty() && (l.back() == '\r' || l.back() == ' ' || l.back() == '\t')) l.remove_suffix(1);
      while (!l.empty() && (l.front() == ' ' || l.front() == '\t')) l.remove_prefix(1);
      if (!l.empty() && l.front() != '#') lines_.emplace_back(l, line_no);
      p = e + 1;
    }
  }

  Model parse();

 private:
  [[noreturn]] void error(const std::string& msg) const {
    const size_t line = k_ > 0 && k_ <= lines_.size() ? lines_[k_ - 1].second : 0;
    throw CbfError(source_ + ":" + std::to_string(line) + ": " + msg);
  }
  // The next data line split into fields.
  std::vector<std::string_view> fields(size_t expect_min, const char* what) {
    if (k_ >= lines_.size()) error(std::string("unexpected end of file in ") + what);
    std::string_view l = lines_[k_++].first;
    std::vector<std::string_view> f;
    while (!l.empty()) {
      size_t s = 0;
      while (s < l.size() && (l[s] == ' ' || l[s] == '\t')) ++s;
      size_t e = s;
      while (e < l.size() && l[e] != ' ' && l[e] != '\t') ++e;
      if (e > s) f.push_back(l.substr(s, e - s));
      l.remove_prefix(e);
    }
    if (f.size() < expect_min) error(std::string("too few fields in ") + what);
    return f;
  }
  long long integer(std::string_view s) const {
    const std::string t(s);
    char* end = nullptr;
    errno = 0;
    const long long v = std::strtoll(t.c_str(), &end, 10);
    if (end != t.c_str() + t.size() || errno != 0) error("invalid integer '" + t + "'");
    return v;
  }
  double number(std::string_view s) const {
    const std::string t(s);
    char* end = nullptr;
    const double v = std::strtod(t.c_str(), &end);
    if (end != t.c_str() + t.size() || !std::isfinite(v)) error("invalid number '" + t + "'");
    return v;
  }
  Index index(std::string_view s, Index n, const char* what) const {
    const long long v = integer(s);
    if (v < 0 || v >= n) error(std::string(what) + " index " + std::string(s) + " out of range");
    return static_cast<Index>(v);
  }
  void cones(std::vector<ConeChunk>& out, Index& total, const char* what) {
    const auto head = fields(2, what);
    total = static_cast<Index>(integer(head[0]));
    const long long k = integer(head[1]);
    if (total < 0 || k < 0) error(std::string("invalid ") + what + " header");
    Index sum = 0;
    for (long long c = 0; c < k; ++c) {
      const auto f = fields(2, what);
      ConeChunk ch{std::string(f[0]), static_cast<Index>(integer(f[1]))};
      if (ch.kind != "F" && ch.kind != "L+" && ch.kind != "L-" && ch.kind != "L=" && ch.kind != "Q" && ch.kind != "QR") {
        error("cone '" + ch.kind + "' is not supported (only F, L+, L-, L=, Q, QR)");
      }
      if (ch.size < 1 || (ch.kind == "QR" && ch.size < 2)) error("invalid cone size");
      sum += ch.size;
      out.push_back(std::move(ch));
    }
    if (sum != total) error(std::string(what) + ": the cones cover " + std::to_string(sum) + " of " + std::to_string(total));
  }

  std::string source_;
  std::vector<std::pair<std::string_view, size_t>> lines_;
  size_t k_ = 0;
};

Model CbfParser::parse() {
  ObjSense sense = ObjSense::Minimize;
  Index n = -1, m = 0;
  std::vector<ConeChunk> var_cones, con_cones;
  std::vector<Index> ints;
  std::vector<std::pair<Index, double>> obj_a;
  double obj_b = 0.0;
  std::vector<std::tuple<Index, Index, double>> a_coord;
  std::vector<std::pair<Index, double>> b_coord;
  while (k_ < lines_.size()) {
    const std::string kw(lines_[k_++].first);
    if (kw == "VER") {
      const long long v = integer(fields(1, "VER")[0]);
      if (v < 1 || v > 3) error("CBF version " + std::to_string(v) + " is not supported (1-3)");
    } else if (kw == "OBJSENSE") {
      const auto f = fields(1, "OBJSENSE");
      if (f[0] == "MIN") sense = ObjSense::Minimize;
      else if (f[0] == "MAX") sense = ObjSense::Maximize;
      else error("OBJSENSE must be MIN or MAX");
    } else if (kw == "VAR") {
      cones(var_cones, n, "VAR");
    } else if (kw == "CON") {
      cones(con_cones, m, "CON");
    } else if (kw == "INT") {
      if (n < 0) error("INT before VAR");
      const long long c = integer(fields(1, "INT")[0]);
      for (long long t = 0; t < c; ++t) ints.push_back(index(fields(1, "INT")[0], n, "INT"));
    } else if (kw == "OBJACOORD") {
      if (n < 0) error("OBJACOORD before VAR");
      const long long c = integer(fields(1, "OBJACOORD")[0]);
      for (long long t = 0; t < c; ++t) {
        const auto f = fields(2, "OBJACOORD");
        obj_a.emplace_back(index(f[0], n, "OBJACOORD"), number(f[1]));
      }
    } else if (kw == "OBJBCOORD") {
      obj_b = number(fields(1, "OBJBCOORD")[0]);
    } else if (kw == "ACOORD") {
      if (n < 0) error("ACOORD before VAR");
      const long long c = integer(fields(1, "ACOORD")[0]);
      for (long long t = 0; t < c; ++t) {
        const auto f = fields(3, "ACOORD");
        a_coord.emplace_back(index(f[0], m, "ACOORD row"), index(f[1], n, "ACOORD column"), number(f[2]));
      }
    } else if (kw == "BCOORD") {
      const long long c = integer(fields(1, "BCOORD")[0]);
      for (long long t = 0; t < c; ++t) {
        const auto f = fields(2, "BCOORD");
        b_coord.emplace_back(index(f[0], m, "BCOORD"), number(f[1]));
      }
    } else {
      error("section '" + kw + "' is not supported (semidefinite and power / exponential cones are not)");
    }
  }
  if (n < 0) error("no VAR section");

  // Columns: bounds from the variable cones.
  std::vector<double> lb(static_cast<size_t>(n), -kInf), ub(static_cast<size_t>(n), kInf);
  std::vector<std::pair<std::string, std::vector<Index>>> cone_members;  // variable and row cones, as columns
  Index j0 = 0;
  for (const ConeChunk& c : var_cones) {
    std::vector<Index> idx;
    for (Index t = 0; t < c.size; ++t) idx.push_back(j0 + t);
    j0 += c.size;
    if (c.kind == "L+") {
      for (Index j : idx) lb[j] = 0.0;
    } else if (c.kind == "L-") {
      for (Index j : idx) ub[j] = 0.0;
    } else if (c.kind == "L=") {
      for (Index j : idx) lb[j] = ub[j] = 0.0;
    } else if (c.kind == "Q" || c.kind == "QR") {
      for (Index t = 0; t < (c.kind == "Q" ? 1 : 2); ++t) lb[idx[t]] = std::max(lb[idx[t]], 0.0);
      cone_members.emplace_back(c.kind, std::move(idx));
    }
  }
  // Rows of A x + b (merged duplicates).
  std::vector<std::vector<std::pair<Index, double>>> rows(static_cast<size_t>(m));
  for (const auto& [i, j, a] : a_coord) rows[i].emplace_back(j, a);
  for (auto& r : rows) {
    std::sort(r.begin(), r.end());
    size_t w = 0;
    for (size_t t = 0; t < r.size(); ++t) {
      if (w > 0 && r[w - 1].first == r[t].first) r[w - 1].second += r[t].second;
      else r[w++] = r[t];
    }
    r.resize(w);
    r.erase(std::remove_if(r.begin(), r.end(), [](const auto& e) { return e.second == 0.0; }), r.end());
  }
  std::vector<double> b(static_cast<size_t>(m), 0.0);
  for (const auto& [i, v] : b_coord) b[i] += v;

  ModelBuilder mb;
  std::vector<char> is_int(static_cast<size_t>(n), 0);
  for (Index j : ints) is_int[j] = 1;
  // Linear rows first, their single-column ones as bounds (before the columns are added).
  struct LinRow {
    Index i;
    std::string kind;
  };
  std::vector<LinRow> lin;
  struct RowCone {
    std::string kind;
    std::vector<Index> rows;
  };
  std::vector<RowCone> row_cones;
  Index i0 = 0;
  for (const ConeChunk& c : con_cones) {
    if (c.kind == "Q" || c.kind == "QR") {
      RowCone rc{c.kind, {}};
      for (Index t = 0; t < c.size; ++t) rc.rows.push_back(i0 + t);
      row_cones.push_back(std::move(rc));
    } else if (c.kind != "F") {
      for (Index t = 0; t < c.size; ++t) {
        const Index i = i0 + t;
        if (rows[i].size() == 1) {
          // a x_j + b (cone) 0: a bound on x_j (flipped for a < 0)
          const auto [j, a] = rows[i][0];
          const double v = -b[i] / a;
          if (c.kind == "L=") {
            lb[j] = std::max(lb[j], v), ub[j] = std::min(ub[j], v);
          } else if ((c.kind == "L+") == (a > 0)) {
            lb[j] = std::max(lb[j], v);
          } else {
            ub[j] = std::min(ub[j], v);
          }
        } else {
          lin.push_back({i, c.kind});
        }
      }
    }
    i0 += c.size;
  }
  std::vector<double> obj(static_cast<size_t>(n), 0.0);
  for (const auto& [j, a] : obj_a) obj[j] += a;
  for (Index j = 0; j < n; ++j) {
    double lo = lb[j], hi = ub[j];
    if (is_int[j]) {
      if (is_finite(lo)) lo = std::ceil(lo - 1e-9);
      if (is_finite(hi)) hi = std::floor(hi + 1e-9);
    }
    mb.add_col("x" + std::to_string(j), obj[j], lo, hi, is_int[j] ? VarType::Integer : VarType::Continuous);
  }
  for (const LinRow& r : lin) {
    // a x + b (cone) 0: L+ a x >= -b, L- a x <= -b, L= a x = -b.
    const double rhs = -b[r.i];
    const Index row = mb.add_row("r" + std::to_string(r.i), r.kind == "L-" ? -kInf : rhs, r.kind == "L+" ? kInf : rhs);
    for (const auto& [j, a] : rows[r.i]) mb.add_coef(row, j, a);
  }
  // Row cones: one column per member (its own column, or a new one defined by an equality row).
  for (const RowCone& rc : row_cones) {
    std::vector<Index> members;
    const Index heads = rc.kind == "Q" ? 1 : 2;
    for (size_t pos = 0; pos < rc.rows.size(); ++pos) {
      const Index i = rc.rows[pos];
      if (b[i] == 0.0 && rows[i].size() == 1 && rows[i][0].second == 1.0) {
        const Index j = rows[i][0].first;
        if (static_cast<Index>(pos) < heads) mb.model().col_lb[j] = std::max(mb.model().col_lb[j], 0.0);
        members.push_back(j);
      } else {
        const Index y = mb.add_col("y" + std::to_string(i), 0.0, static_cast<Index>(pos) < heads ? 0.0 : -kInf, kInf);
        const Index d = mb.add_row("d" + std::to_string(i), b[i], b[i]);  // y - a x = b
        mb.add_coef(d, y, 1.0);
        for (const auto& [j, a] : rows[i]) mb.add_coef(d, j, -a);
        members.push_back(y);
      }
    }
    cone_members.emplace_back(rc.kind, std::move(members));
  }
  Model model = mb.build();
  model.sense = sense;
  model.obj_offset = obj_b;
  Index kcone = 0;
  for (const auto& [kind, idx] : cone_members) {
    const size_t heads = kind == "Q" ? 1 : 2;
    if (idx.size() <= heads) continue;  // only the heads' signs (bounds)
    model.row_names.push_back("k" + std::to_string(kcone++));
    model.row_lb.push_back(-kInf);
    model.row_ub.push_back(0.0);
    model.A.num_minor = model.num_rows();
    const Index row = model.num_rows() - 1;
    for (size_t t = heads; t < idx.size(); ++t) model.add_qrow(row, idx[t], idx[t], 1.0);
    if (kind == "Q") model.add_qrow(row, idx[0], idx[0], -1.0);
    else model.add_qrow(row, idx[0], idx[1], -2.0);
  }
  model.compact_qrows();
  return model;
}

}  // namespace

Model read_cbf_string(const std::string& text, const std::string& source) {
  CbfParser parser(text, source);
  return parser.parse();
}

Model read_cbf(const std::string& path) {
  std::string text;
  try {
    text = read_maybe_gzipped_file(path);
  } catch (const std::bad_alloc&) {
    throw;
  } catch (const std::exception& e) {
    throw CbfError(e.what());
  }
  return read_cbf_string(text, path);
}

}  // namespace milp
