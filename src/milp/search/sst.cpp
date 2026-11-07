#include "milp/search/sst.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace milp {

namespace {

Index uf_find(std::vector<Index>& uf, Index j) {
  while (uf[j] != j) j = uf[j] = uf[uf[j]];
  return j;
}

// Row activity bounds over the column bounds (finite parts and counts of infinite contributions).
struct RowActivity {
  std::vector<double> min_fin, max_fin;
  std::vector<int> min_inf, max_inf;
};

double min_contrib(double a, double lb, double ub) { return a > 0 ? a * lb : a * ub; }
double max_contrib(double a, double lb, double ub) { return a > 0 ? a * ub : a * lb; }

RowActivity row_activity(const MipProblem& p) {
  const Index m = p.num_rows();
  RowActivity act{std::vector<double>(m, 0.0), std::vector<double>(m, 0.0), std::vector<int>(m, 0),
                  std::vector<int>(m, 0)};
  for (Index j = 0; j < p.num_cols(); ++j) {
    for (Index k = p.A.begin(j); k < p.A.end(j); ++k) {
      const Index r = p.A.index[k];
      const double a = p.A.value[k];
      const double lo = min_contrib(a, p.col_lb[j], p.col_ub[j]), hi = max_contrib(a, p.col_lb[j], p.col_ub[j]);
      if (is_finite(lo)) act.min_fin[r] += lo;
      else ++act.min_inf[r];
      if (is_finite(hi)) act.max_fin[r] += hi;
      else ++act.max_inf[r];
    }
  }
  return act;
}

uint64_t mix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

// Generators of a subgroup of the stabilizer of l in <gens>: the generators fixing l plus random
// Schreier generators u_z^-1 s u_x (exact products of group elements, each fixing l). Works on
// the support of the generators only. Deterministic for a given seed; adds the steps to work.
std::vector<std::vector<std::pair<Index, Index>>> stabilizer(const std::vector<std::vector<std::pair<Index, Index>>>& gens,
                                                             Index l, Index n, uint64_t seed, long long& work,
                                                             int max_kept, int max_new) {
  std::vector<std::vector<std::pair<Index, Index>>> out;
  std::vector<int> tg;
  auto image = [&](int g, Index j) {
    const auto& v = gens[g];
    auto it = std::lower_bound(v.begin(), v.end(), j, [](const std::pair<Index, Index>& e, Index c) { return e.first < c; });
    return it != v.end() && it->first == j ? it->second : j;
  };
  for (int g = 0; g < static_cast<int>(gens.size()); ++g) {
    work += 8;
    if (image(g, l) == l) {
      if (static_cast<int>(out.size()) < max_kept) out.push_back(gens[g]);
    } else if (tg.size() < 32) {
      tg.push_back(g);
    }
  }
  if (tg.empty()) return out;
  // Local numbering of the support of the tree generators.
  std::vector<Index> supp;
  for (int g : tg) {
    for (const auto& e : gens[g]) supp.push_back(e.first);
  }
  std::sort(supp.begin(), supp.end());
  supp.erase(std::unique(supp.begin(), supp.end()), supp.end());
  const Index s = static_cast<Index>(supp.size());
  auto loc = [&](Index j) { return static_cast<Index>(std::lower_bound(supp.begin(), supp.end(), j) - supp.begin()); };
  const int T = static_cast<int>(tg.size());
  std::vector<std::vector<Index>> d(T, std::vector<Index>(s));
  for (int k = 0; k < T; ++k) {
    for (Index i = 0; i < s; ++i) d[k][i] = i;
    for (const auto& [a, b] : gens[tg[k]]) d[k][loc(a)] = loc(b);
    work += s + 4 * static_cast<long long>(gens[tg[k]].size());
  }
  const Index L = loc(l);
  std::vector<Index> parent(s, -1), orbit(1, L);
  std::vector<int> via(s, -1), depth(s, 0);
  parent[L] = L;
  for (size_t h = 0; h < orbit.size(); ++h) {
    const Index x = orbit[h];
    for (int k = 0; k < T; ++k) {
      const Index z = d[k][x];
      if (parent[z] >= 0) continue;
      parent[z] = x;
      via[z] = k;
      depth[z] = depth[x] + 1;
      orbit.push_back(z);
    }
  }
  work += static_cast<long long>(orbit.size()) * T;
  std::vector<Index> ux, uz, uzinv(s);
  std::vector<int> path;
  auto word = [&](Index x, std::vector<Index>& u) {
    u.resize(s);
    for (Index c = 0; c < s; ++c) u[c] = c;
    path.clear();
    for (; x != L; x = parent[x]) path.push_back(via[x]);
    for (auto it = path.rbegin(); it != path.rend(); ++it) {
      const std::vector<Index>& p = d[*it];
      for (Index c = 0; c < s; ++c) u[c] = p[u[c]];
    }
    work += static_cast<long long>(s) * (static_cast<long long>(path.size()) + 1);
  };
  uint64_t rng = mix64(seed);
  int added = 0;
  for (int a = 0; a < 3 * max_new && added < max_new; ++a) {
    rng = mix64(rng);
    const Index x = orbit[rng % orbit.size()];
    rng = mix64(rng);
    const int k = static_cast<int>(rng % static_cast<uint64_t>(T));
    const Index z = d[k][x];
    word(x, ux);
    word(z, uz);
    for (Index c = 0; c < s; ++c) uzinv[uz[c]] = c;
    std::vector<std::pair<Index, Index>> g;
    for (Index c = 0; c < s; ++c) {
      const Index v = uzinv[d[k][ux[c]]];
      if (v != c) g.emplace_back(supp[c], supp[v]);
    }
    work += 2 * static_cast<long long>(s);
    if (g.empty() || std::find(out.begin(), out.end(), g) != out.end()) continue;
    out.push_back(std::move(g));
    ++added;
  }
  (void)n;
  return out;
}

class Sst {
 public:
  Sst(const MipProblem& p, const SstOptions& o) : p_(p), o_(o), n_(p.num_cols()) {}
  SstResult run(const Symmetry& sym, const std::vector<char>& skip_gen);

 private:
  bool binary(Index j) const { return p_.is_int[j] && p_.col_lb[j] == 0.0 && p_.col_ub[j] == 1.0; }
  // Marks (in conf_) the binary members of the marked orbit that cannot be 1 together with l.
  int conflicts(Index l);
  bool pair_conflict(Index r, double al, Index l, double af, Index f) const;

  const MipProblem& p_;
  SstOptions o_;
  Index n_;
  SparseMatrix R_;  // CSR
  RowActivity act_;
  std::vector<char> in_orbit_, conf_;
  std::vector<Index> conf_list_;
  long long work_ = 0;
};

bool Sst::pair_conflict(Index r, double al, Index l, double af, Index f) const {
  // Activity range of row r with x_l = x_f = 1 (both binary).
  const double lbl = p_.col_lb[l], ubl = p_.col_ub[l], lbf = p_.col_lb[f], ubf = p_.col_ub[f];
  const double ub = p_.row_ub[r], lb = p_.row_lb[r];
  if (is_finite(ub) && act_.min_inf[r] == 0) {
    const double mn = act_.min_fin[r] - min_contrib(al, lbl, ubl) - min_contrib(af, lbf, ubf) + al + af;
    if (mn > ub + 1e-6 + 1e-9 * std::abs(ub)) return true;
  }
  if (is_finite(lb) && act_.max_inf[r] == 0) {
    const double mx = act_.max_fin[r] - max_contrib(al, lbl, ubl) - max_contrib(af, lbf, ubf) + al + af;
    if (mx < lb - 1e-6 - 1e-9 * std::abs(lb)) return true;
  }
  return false;
}

int Sst::conflicts(Index l) {
  int cnt = 0;
  if (!binary(l)) return 0;
  for (Index k = p_.A.begin(l); k < p_.A.end(l); ++k) {
    const Index r = p_.A.index[k];
    const double al = p_.A.value[k];
    work_ += R_.end(r) - R_.begin(r) + 1;
    for (Index e = R_.begin(r); e < R_.end(r); ++e) {
      const Index f = R_.index[e];
      if (f == l || !in_orbit_[f] || conf_[f] || !binary(f)) continue;
      if (pair_conflict(r, al, l, R_.value[e], f)) {
        conf_[f] = 1;
        conf_list_.push_back(f);
        ++cnt;
      }
    }
  }
  return cnt;
}

SstResult Sst::run(const Symmetry& sym, const std::vector<char>& skip_gen) {
  SstResult out;
  std::vector<std::vector<std::pair<Index, Index>>> H;
  for (size_t g = 0; g < sym.gens.size(); ++g) {
    if (skip_gen.empty() || !skip_gen[g]) H.push_back(sym.gens[g]);
  }
  if (H.empty() || n_ == 0) return out;
  R_ = p_.A.transpose();
  act_ = row_activity(p_);
  const long long size = static_cast<long long>(n_) + p_.num_rows() + p_.A.nnz();
  const long long work_limit = o_.work_limit > 0 ? o_.work_limit : std::min<long long>(40 * size + 20000000LL, 300000000LL);
  const long long max_rows = o_.max_rows >= 0 ? o_.max_rows : 2LL * n_ + 1000;
  work_ = 2 * size;
  in_orbit_.assign(n_, 0);
  conf_.assign(n_, 0);
  std::vector<char> zero(n_, 0), leader(n_, 0);
  std::vector<Index> uf(n_), members;
  std::vector<Index> head(n_, -1), next(n_, -1);  // orbit member lists by root
  std::vector<int> score(n_, -1);
  std::vector<Index> score_size(n_, 0);
  for (int step = 0; step < o_.max_leaders && !H.empty(); ++step) {
    // Orbits of H.
    std::iota(uf.begin(), uf.end(), 0);
    for (const auto& g : H) {
      work_ += static_cast<long long>(g.size());
      for (const auto& [a, b] : g) {
        const Index ra = uf_find(uf, a), rb = uf_find(uf, b);
        if (ra != rb) uf[std::max(ra, rb)] = std::min(ra, rb);
      }
    }
    work_ += n_;
    for (Index j = n_; j-- > 0;) {
      const Index r = uf_find(uf, j);
      if (r == j && head[j] < 0) continue;
      next[j] = head[r];
      head[r] = j;
    }
    // Members are listed in increasing order from head[root] (root = smallest member).
    Index best = -1;
    int best_conf = -1;
    size_t best_size = 0;
    bool fixed_more = false;
    for (Index r = 0; r < n_; ++r) {
      if (work_ > work_limit) break;  // scoring the orbits alone can take long (bab2: long rows, minutes)
      if (uf[r] != r || head[r] < 0) continue;
      // Symmetries map columns onto columns of the same type and bounds: the root stands for its orbit.
      if (o_.mode == 3 ? binary(r) : !p_.is_int[r]) continue;
      members.clear();
      bool has_zero = false;
      for (Index j = r; j >= 0; j = next[j]) {
        members.push_back(j);
        has_zero |= zero[j] != 0;
      }
      if (members.size() < 2) continue;
      if (has_zero) {
        // The orbit of a column fixed to 0 is 0 (the reduced problem is H-invariant).
        for (Index j : members) {
          if (!zero[j] && binary(j)) {
            zero[j] = 1;
            out.fixed_zero.push_back(j);
            fixed_more = true;
          }
        }
        continue;
      }
      // Orbits of the stabilizer refine the previous ones: an orbit with the same smallest member
      // and size is the same set, with the same conflict count.
      int c = o_.mode == 3 ? 0 : score_size[r] == static_cast<Index>(members.size()) ? score[r] : -1;
      if (c < 0) {
        for (Index j : members) in_orbit_[j] = 1;
        c = conflicts(r);
        for (Index j : members) in_orbit_[j] = 0;
        for (Index f : conf_list_) conf_[f] = 0;
        conf_list_.clear();
        score[r] = c;
        score_size[r] = static_cast<Index>(members.size());
      }
      if (c > best_conf || (c == best_conf && members.size() > best_size)) {
        best = r;
        best_conf = c;
        best_size = members.size();
      }
    }
    if (work_ > work_limit) {
      out.budget = true;
      break;
    }
    if (best < 0 || (o_.mode == 1 && best_conf <= 0)) {
      for (Index j = 0; j < n_; ++j) head[j] = next[j] = -1;
      if (fixed_more) continue;  // orbits fixed to 0 do not change H; nothing else to do
      break;
    }
    // Leader: the smallest column of the chosen orbit.
    const Index l = best;
    members.clear();
    for (Index j = l; j >= 0; j = next[j]) members.push_back(j);
    for (Index j = 0; j < n_; ++j) head[j] = next[j] = -1;
    for (Index j : members) in_orbit_[j] = 1;
    conflicts(l);
    for (Index j : members) in_orbit_[j] = 0;
    for (Index f : members) {
      if (f == l) continue;
      if (conf_[f]) {
        zero[f] = 1;
        out.fixed_zero.push_back(f);
      } else {
        out.rows.emplace_back(l, f);
      }
    }
    for (Index f : conf_list_) conf_[f] = 0;
    conf_list_.clear();
    leader[l] = 1;
    ++out.leaders;
    if (static_cast<long long>(out.rows.size()) > max_rows || work_ > work_limit) {
      out.budget = true;
      break;
    }
    // H := (a subgroup of) the stabilizer of l in H.
    const uint64_t seed = o_.seed * 0x9e3779b97f4a7c15ULL + static_cast<uint64_t>(step) * 0xbf58476d1ce4e5b9ULL + 1;
    H = stabilizer(H, l, n_, seed, work_, 64, 16);
    if (work_ > work_limit) {
      out.budget = true;
      break;
    }
  }
  std::sort(out.fixed_zero.begin(), out.fixed_zero.end());
  out.fixed_zero.erase(std::unique(out.fixed_zero.begin(), out.fixed_zero.end()), out.fixed_zero.end());
  // Rows whose follower was fixed to 0 later are redundant.
  out.rows.erase(std::remove_if(out.rows.begin(), out.rows.end(), [&](const std::pair<Index, Index>& e) { return zero[e.second] != 0; }),
                 out.rows.end());
  // Without rows, the fixings alone are a relaxation of the full system (they are implied by
  // it), so an optimum still survives.
  if (!o_.rows) out.rows.clear();
  out.work = work_;
  return out;
}

}  // namespace

SstResult sst_reductions(const MipProblem& p, const Symmetry& sym, const std::vector<char>& skip_gen,
                         const SstOptions& opts) {
  Sst s(p, opts);
  return s.run(sym, skip_gen);
}

MipProblem apply_sst(const MipProblem& p, const SstResult& r) {
  MipProblem q = p;
  for (Index j : r.fixed_zero) q.col_ub[j] = 0.0;
  if (r.rows.empty()) return q;
  const Index m0 = p.num_rows(), n = p.num_cols();
  std::vector<Index> J, I;
  std::vector<double> V;
  const size_t nnz = static_cast<size_t>(p.A.nnz()) + 2 * r.rows.size();
  J.reserve(nnz);
  I.reserve(nnz);
  V.reserve(nnz);
  for (Index j = 0; j < n; ++j) {
    for (Index k = p.A.begin(j); k < p.A.end(j); ++k) J.push_back(j), I.push_back(p.A.index[k]), V.push_back(p.A.value[k]);
  }
  for (size_t t = 0; t < r.rows.size(); ++t) {
    const Index row = m0 + static_cast<Index>(t);
    J.push_back(r.rows[t].first), I.push_back(row), V.push_back(1.0);
    J.push_back(r.rows[t].second), I.push_back(row), V.push_back(-1.0);
  }
  q.A = SparseMatrix::from_triplets(n, m0 + static_cast<Index>(r.rows.size()), J, I, V);
  q.row_lb.resize(m0 + r.rows.size(), 0.0);
  q.row_ub.resize(m0 + r.rows.size(), kInf);
  return q;
}

}  // namespace milp
