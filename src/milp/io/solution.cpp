#include "milp/io/solution.hpp"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

#include "milp/io/charconv_compat.hpp"

namespace milp {

std::string format_solution(const Model& model, const SolveResult& result) {
  std::string out;
  out.reserve(64 + (result.has_solution ? 32 * static_cast<size_t>(model.num_cols()) : 0));
  out += "# Solution for model " + model.name + "\n";
  out += std::string("# Status = ") + to_string(result.status) + "\n";
  if (!result.has_solution) return out;
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.17g", result.objective);
  out += "# Objective value = ";
  out += buf;
  out += '\n';
  const size_t n = std::min(result.x.size(), model.col_names.size());
  for (size_t j = 0; j < n; ++j) {
    out += model.col_names[j];
    out += ' ';
    const double v = result.x[j] == 0.0 ? 0.0 : result.x[j];  // no "-0"
    out.append(buf, static_cast<size_t>(io_detail::format_double(buf, v)));
    out += '\n';
  }
  return out;
}

std::string format_duals(const Model& model, const SolveResult& result) {
  std::string out = "# Duals for model " + model.name + "\n";
  out += std::string("# Status = ") + to_string(result.status) + "\n";
  if (!result.has_duals) {
    out += "# No dual values (LPs and continuous QPs solved to optimality, MILPs with a solution)\n";
    return out;
  }
  if (result.fixed_mip_duals) out += "# Fixed-MIP duals: the LP with the integer variables fixed at the solution\n";
  char buf[64];
  auto line = [&](const char* kind, const std::string& name, double v) {
    out += kind;
    out += name;
    out += ' ';
    if (v == 0.0) v = 0.0;  // no "-0"
    out.append(buf, static_cast<size_t>(io_detail::format_double(buf, v)));
    out += '\n';
  };
  for (size_t i = 0; i < result.row_dual.size() && i < model.row_names.size(); ++i) line("row ", model.row_names[i], result.row_dual[i]);
  for (size_t j = 0; j < result.reduced_cost.size() && j < model.col_names.size(); ++j)
    line("col ", model.col_names[j], result.reduced_cost[j]);
  return out;
}

namespace {

void write_text(const std::string& path, const std::string& text) {
  std::FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) throw std::runtime_error(path + ": cannot open file for writing");
  const bool wrote = std::fwrite(text.data(), 1, text.size(), f) == text.size();
  const bool closed = std::fclose(f) == 0;
  if (!wrote || !closed) throw std::runtime_error(path + ": write error");
}

}  // namespace

void write_duals(const std::string& path, const Model& model, const SolveResult& result) {
  write_text(path, format_duals(model, result));
}

void write_solution(const std::string& path, const Model& model, const SolveResult& result) {
  write_text(path, format_solution(model, result));
}

std::vector<double> read_solution_string(const std::string& text, const Model& model) {
  std::unordered_map<std::string_view, Index> index;
  index.reserve(model.col_names.size());
  for (size_t j = 0; j < model.col_names.size(); ++j) index.emplace(model.col_names[j], static_cast<Index>(j));

  std::vector<double> x(model.col_names.size(), std::numeric_limits<double>::quiet_NaN());
  auto space = [](char c) { return c == ' ' || c == '\t' || c == '\r'; };
  const char* p = text.data();
  const char* end = p + text.size();
  for (size_t line = 1; p < end; ++line) {
    const char* nl = static_cast<const char*>(std::memchr(p, '\n', static_cast<size_t>(end - p)));
    const char* le = nl ? nl : end;
    std::string_view tok[3];
    int ntok = 0;
    for (const char* q = p; q < le;) {
      while (q < le && space(*q)) ++q;
      if (q == le) break;
      const char* s = q;
      while (q < le && !space(*q)) ++q;
      if (ntok < 3) tok[ntok] = std::string_view(s, static_cast<size_t>(q - s));
      ++ntok;
    }
    p = nl ? nl + 1 : end;
    if (ntok == 0 || tok[0].front() == '#') continue;
    const std::string where = "solution line " + std::to_string(line) + ": ";
    if (ntok != 2) throw std::runtime_error(where + "expected '<name> <value>'");
    auto it = index.find(tok[0]);
    if (it == index.end()) throw std::runtime_error(where + "unknown variable '" + std::string(tok[0]) + "'");
    double v;
    if (!io_detail::parse_double(tok[1], v)) throw std::runtime_error(where + "invalid number '" + std::string(tok[1]) + "'");
    x[it->second] = v;
  }
  return x;
}

std::vector<double> read_solution(const std::string& path, const Model& model) {
  std::FILE* f = std::fopen(path.c_str(), "rb");
  if (!f) throw std::runtime_error(path + ": cannot open file");
  std::string data;
  char buf[1 << 16];
  size_t got;
  while ((got = std::fread(buf, 1, sizeof(buf), f)) > 0) data.append(buf, got);
  std::fclose(f);
  try {
    return read_solution_string(data, model);
  } catch (const std::exception& e) {
    throw std::runtime_error(path + ": " + e.what());
  }
}

}  // namespace milp
