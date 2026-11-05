#include "milp/search/orbitope.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <numeric>
#include <unordered_map>

namespace milp {

namespace {

uint64_t mix64(uint64_t x) {
  x += 0x9e3779b97f4a7c15ULL;
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

uint64_t hash_double(double v) {
  if (v == 0.0) v = 0.0;  // -0 == +0
  uint64_t u;
  std::memcpy(&u, &v, sizeof u);
  return mix64(u);
}

Index find_root(std::vector<Index>& parent, Index a) {
  while (parent[a] != a) a = parent[a] = parent[parent[a]];
  return a;
}

// Exact check that a column permutation (with some row permutation) maps the problem onto itself.
class PermChecker {
 public:
  PermChecker(const MipProblem& p, const SparseMatrix& R) : p_(p), R_(R) {
    const Index m = p.num_rows();
    hash_.resize(m);
    used_.assign(m, 0);
    for (Index r = 0; r < m; ++r) {
      hash_[r] = row_hash(r, nullptr);
      by_hash_.emplace(hash_[r], r);
    }
  }
  // sigma: dense column permutation (identity outside `moved`).
  bool check(const std::vector<Index>& sigma, const std::vector<Index>& moved) {
    for (Index j : moved) {
      const Index k = sigma[j];
      if (p_.cost[j] != p_.cost[k] || p_.col_lb[j] != p_.col_lb[k] || p_.col_ub[j] != p_.col_ub[k] ||
          p_.is_int[j] != p_.is_int[k]) {
        return false;
      }
    }
    // Every row touching a moved column must map onto a distinct row with the same data.
    std::vector<Index> rows;
    for (Index j : moved) {
      for (Index k = p_.A.begin(j); k < p_.A.end(j); ++k) rows.push_back(p_.A.index[k]);
    }
    std::sort(rows.begin(), rows.end());
    rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
    // Each row must map onto a distinct row with the same content (duplicate rows map onto
    // distinct copies, so the row multiset is invariant).
    ++stamp_;
    for (Index r : rows) {
      const uint64_t h = row_hash(r, &sigma);
      Index img = -1;
      auto range = by_hash_.equal_range(h);
      for (auto it = range.first; it != range.second && img < 0; ++it) {
        if (used_[it->second] != stamp_ && same_row(r, it->second, sigma)) {
          img = it->second;
        }
      }
      if (img < 0) return false;
      used_[img] = stamp_;
    }
    return true;
  }

 private:
  uint64_t row_hash(Index r, const std::vector<Index>* sigma) const {
    uint64_t h = mix64(hash_double(p_.row_lb[r]) ^ (hash_double(p_.row_ub[r]) * 31));
    for (Index k = R_.begin(r); k < R_.end(r); ++k) {
      const Index j = sigma ? (*sigma)[R_.index[k]] : R_.index[k];
      h += mix64(static_cast<uint64_t>(j) * 0x9e3779b97f4a7c15ULL ^ hash_double(R_.value[k]));
    }
    return h;
  }
  bool same_row(Index r, Index t, const std::vector<Index>& sigma) {
    if (p_.row_lb[r] != p_.row_lb[t] || p_.row_ub[r] != p_.row_ub[t]) return false;
    if (R_.end(r) - R_.begin(r) != R_.end(t) - R_.begin(t)) return false;
    a_.clear();
    b_.clear();
    for (Index k = R_.begin(r); k < R_.end(r); ++k) a_.emplace_back(sigma[R_.index[k]], R_.value[k]);
    for (Index k = R_.begin(t); k < R_.end(t); ++k) b_.emplace_back(R_.index[k], R_.value[k]);
    std::sort(a_.begin(), a_.end());
    std::sort(b_.begin(), b_.end());
    return a_ == b_;
  }
  const MipProblem& p_;
  const SparseMatrix& R_;
  std::vector<uint64_t> hash_;
  std::unordered_multimap<uint64_t, Index> by_hash_;
  std::vector<std::pair<Index, double>> a_, b_;
  std::vector<long long> used_;  // stamp per row: already an image in the current check
  long long stamp_ = 0;
};

// Whether the entries cols (binary) lie in one model row allowing at most one of them to be 1.
bool is_packing(const MipProblem& p, const SparseMatrix& R, const std::vector<Index>& cols, std::vector<char>& mark) {
  for (Index j : cols) {
    if (!p.is_int[j] || p.col_lb[j] < 0.0 || p.col_ub[j] > 1.0) return false;
  }
  for (Index j : cols) mark[j] = 1;
  bool found = false;
  const Index j0 = cols[0];
  for (Index k = p.A.begin(j0); k < p.A.end(j0) && !found; ++k) {
    const Index r = p.A.index[k];
    // A sum of nonnegative terms <= row_ub where any two of the entries already exceed row_ub.
    if (!is_finite(p.row_ub[r])) continue;
    double amin = kInf;
    size_t hit = 0;
    bool nonneg = true;
    for (Index q = R.begin(r); q < R.end(r) && nonneg; ++q) {
      const Index j = R.index[q];
      const double a = R.value[q];
      if (a < 0.0 || p.col_lb[j] < 0.0) nonneg = false;
      if (mark[j]) amin = std::min(amin, a), ++hit;
    }
    found = nonneg && hit == cols.size() && amin > 0.0 &&
            2.0 * amin > p.row_ub[r] + 1e-9 * (1.0 + std::abs(p.row_ub[r]));
  }
  for (Index j : cols) mark[j] = 0;
  return found;
}

}  // namespace

std::vector<Orbitope> detect_packing_orbitopes(const MipProblem& p, const Symmetry& sym,
                                               std::vector<char>* used_gen) {
  std::vector<Orbitope> out;
  const Index n = p.num_cols();
  if (used_gen) used_gen->assign(sym.gens.size(), 0);
  if (sym.gens.empty()) return out;
  static const bool dbg = std::getenv("MILP_ORBITOPE_DEBUG") != nullptr;
  // Components (columns moved together by some generator) and orbits.
  std::vector<Index> comp(n), orb(n);
  std::iota(comp.begin(), comp.end(), 0);
  std::iota(orb.begin(), orb.end(), 0);
  std::vector<char> moved(n, 0);
  for (const auto& gen : sym.gens) {
    for (const auto& [a, b] : gen) {
      moved[a] = moved[b] = 1;
      Index ra = find_root(comp, a), rb = find_root(comp, gen[0].first);  // the whole support
      if (ra != rb) comp[ra] = rb;
      ra = find_root(orb, a), rb = find_root(orb, b);
      if (ra != rb) orb[ra] = rb;
    }
  }
  std::vector<std::vector<Index>> comp_cols(n);
  for (Index j = 0; j < n; ++j) {
    if (moved[j]) comp_cols[find_root(comp, j)].push_back(j);
  }
  std::vector<std::vector<int>> comp_gens(n);
  for (size_t g = 0; g < sym.gens.size(); ++g) {
    if (!sym.gens[g].empty()) comp_gens[find_root(comp, sym.gens[g][0].first)].push_back(static_cast<int>(g));
  }
  const SparseMatrix R = p.A.transpose();
  std::unique_ptr<PermChecker> checker;
  std::vector<Index> orb_id(n, -1), block(n), sigma(n);
  std::iota(sigma.begin(), sigma.end(), 0);
  std::vector<char> mark(n, 0);
  // Work (row entries hashed) for inferring the swaps' images, over all components: per component it
  // let bab2 (many components, 147k columns) spend minutes here, past any time limit.
  long long work = 0;
  static const long long max_work = [] {
    const char* v = std::getenv("MILP_ORBITOPE_WORK");
    return v ? std::atoll(v) : 400000000LL;
  }();
  for (Index root = 0; root < n; ++root) {
    const std::vector<Index>& cols = comp_cols[root];
    if (cols.empty()) continue;
    if (work > max_work) {
      if (dbg) std::fprintf(stderr, "ORB work limit: components left unchecked\n");
      break;
    }
    // Orbits of the component (numbered by their smallest column) must all have the same size p.
    std::vector<std::vector<Index>> orbits;
    for (Index j : cols) {
      const Index o = find_root(orb, j);
      if (orb_id[o] < 0) orb_id[o] = static_cast<Index>(orbits.size()), orbits.emplace_back();
      orbits[orb_id[o]].push_back(j);
    }
    auto cleanup = [&]() {
      for (Index j : cols) orb_id[find_root(orb, j)] = -1;
    };
    // Orbit sizes must be multiples of the smallest one, p (the number of blocks); the orbits of
    // size p ("small") define the blocks, the others are permuted along (e.g. pairs of blocks).
    size_t np = orbits[0].size();
    for (const auto& o : orbits) np = std::min(np, o.size());
    bool ok = np >= 2;
    for (const auto& o : orbits) ok = ok && o.size() % np == 0;
    std::vector<char> small(orbits.size(), 0);
    std::vector<size_t> small_ids;
    for (size_t t = 0; t < orbits.size(); ++t) {
      if (orbits[t].size() == np) small[t] = 1, small_ids.push_back(t);
    }
    if (!ok) {
      if (dbg) std::fprintf(stderr, "ORB comp of %zu cols: orbit sizes not multiples of %zu\n", cols.size(), np);
      cleanup();
      continue;
    }
    auto oid = [&](Index j) { return orb_id[find_root(orb, j)]; };
    auto in_comp = [&](Index j) { return moved[j] && find_root(comp, j) == root; };
    // Blocks (the matrix columns, e.g. one truck): small-orbit columns sharing a model row that
    // meets every small orbit at most once.
    for (Index j : cols) block[j] = j;
    std::vector<Index> crow;
    for (Index j : cols) {
      for (Index k = p.A.begin(j); k < p.A.end(j); ++k) crow.push_back(p.A.index[k]);
    }
    std::sort(crow.begin(), crow.end());
    crow.erase(std::unique(crow.begin(), crow.end()), crow.end());
    for (Index r : crow) {
      Index first = -1;
      bool single = true;
      for (Index k = R.begin(r); k < R.end(r) && single; ++k) {
        const Index j = R.index[k];
        if (!in_comp(j) || !small[oid(j)]) continue;
        if (mark[oid(j)]) single = false;
        mark[oid(j)] = 1;
      }
      for (Index k = R.begin(r); k < R.end(r); ++k) {
        const Index j = R.index[k];
        if (in_comp(j) && small[oid(j)]) mark[oid(j)] = 0;
      }
      if (!single) continue;
      for (Index k = R.begin(r); k < R.end(r); ++k) {
        const Index j = R.index[k];
        if (!in_comp(j) || !small[oid(j)]) continue;
        if (first < 0) {
          first = j;
          continue;
        }
        const Index ra = find_root(block, j), rb = find_root(block, first);
        if (ra != rb) block[ra] = rb;
      }
    }
    // Each block must hold exactly one column of every small orbit.
    std::vector<Index> bid(n, -1);
    std::vector<std::vector<Index>> blocks;  // blocks[b][orbit] = column (small orbits only)
    for (Index j : orbits[small_ids[0]]) {
      const Index b = find_root(block, j);
      if (bid[b] >= 0) {
        ok = false;
        break;
      }
      bid[b] = static_cast<Index>(blocks.size());
      blocks.emplace_back(orbits.size(), -1);
    }
    // Small orbits outside the blocks of the first one (a block split in several parts) are
    // treated like the larger orbits: their images under a swap are inferred.
    {
      std::vector<size_t> keep;
      for (size_t t : small_ids) {
        bool in = true;
        for (Index j : orbits[t]) in = in && bid[find_root(block, j)] >= 0;
        if (in) keep.push_back(t);
        else small[t] = 0;
      }
      small_ids = std::move(keep);
    }
    for (size_t t : small_ids) {
      for (Index j : orbits[t]) {
        if (!ok) break;
        const Index b = bid[find_root(block, j)];
        if (b < 0 || blocks[b][t] >= 0) {
          ok = false;
          break;
        }
        blocks[b][t] = j;
      }
    }
    if (!ok || blocks.size() != np) {
      if (dbg) {
        std::vector<Index> roots;
        for (size_t t : small_ids) for (Index j : orbits[t]) roots.push_back(find_root(block, j));
        std::sort(roots.begin(), roots.end());
        const size_t nb = std::unique(roots.begin(), roots.end()) - roots.begin();
        std::fprintf(stderr, "ORB comp of %zu cols, %zu orbits (%zu small), p %zu: no block structure (%zu blocks, %zu from orbit 0)\n",
                     cols.size(), orbits.size(), small_ids.size(), np, nb, blocks.size());
      }
      cleanup();
      continue;
    }
    // Swapping two blocks permutes the other orbits' columns too (a block transposition is an
    // involution). Their images are inferred by hashing each column's rows: resolved columns by
    // their image, unresolved ones by orbit only (the unresolved set is closed under the swap, so
    // both sides hash alike). Columns with a unique candidate image are resolved first, in
    // pairs; then those that look unchanged stay. Everything is verified exactly afterwards.
    std::vector<Index> big;
    for (size_t t = 0; t < orbits.size(); ++t) {
      if (!small[t]) big.insert(big.end(), orbits[t].begin(), orbits[t].end());
    }
    std::vector<char> unres(n, 0);
    auto key = [&](Index k, double a, bool image) -> uint64_t {
      if (unres[k]) return mix64((0xabcdef12345ULL + static_cast<uint64_t>(oid(k))) ^ hash_double(a));
      const Index fk = image ? sigma[k] : k;
      return mix64(static_cast<uint64_t>(fk) * 0x9e3779b97f4a7c15ULL ^ hash_double(a));
    };
    auto col_hash = [&](Index j, bool image) {
      uint64_t h = 0;
      for (Index e = p.A.begin(j); e < p.A.end(j); ++e) {
        const Index r = p.A.index[e];
        uint64_t hr = mix64(hash_double(p.A.value[e]) ^ hash_double(p.row_lb[r]) * 7 ^ hash_double(p.row_ub[r]) * 13);
        for (Index k = R.begin(r); k < R.end(r); ++k) {
          if (R.index[k] != j) hr += key(R.index[k], R.value[k], image);
        }
        work += R.end(r) - R.begin(r);
        h += mix64(hr);
      }
      return h;
    };
    if (!checker) checker = std::make_unique<PermChecker>(p, R);
    std::vector<Index> mv;
    std::vector<uint64_t> hid(n, 0), hsg(n, 0);
    std::vector<std::vector<std::pair<Index, Index>>> swaps;
    for (size_t b = 0; b + 1 < blocks.size() && ok; ++b) {
      mv.clear();
      for (size_t t : small_ids) {
        const Index u = blocks[b][t], v = blocks[b + 1][t];
        sigma[u] = v, sigma[v] = u;
        mv.push_back(u), mv.push_back(v);
      }
      for (Index j : big) unres[j] = 1;
      size_t left = big.size();
      bool lenient = false;
      while (left > 0 && ok) {
        if (work > max_work) {
          if (dbg) std::fprintf(stderr, "ORB   swap %zu: work limit, %zu left\n", b, left);
          ok = false;
          break;
        }
        std::unordered_multimap<uint64_t, Index> by_h;
        for (Index j : big) {
          if (!unres[j]) continue;
          if (work > max_work) break;  // one pass alone can exceed the budget
          hid[j] = col_hash(j, false);
          hsg[j] = col_hash(j, true);
          by_h.emplace(hid[j], j);
        }
        if (work > max_work) continue;
        size_t progress = 0;
        for (Index j : big) {
          if (!unres[j]) continue;
          if (hsg[j] == hid[j]) continue;
          Index img = -1;
          int cnt = 0;
          auto range = by_h.equal_range(hsg[j]);
          for (auto it = range.first; it != range.second; ++it) {
            if (oid(it->second) == oid(j) && unres[it->second]) img = it->second, ++cnt;
          }
          if (cnt != 1 || img == j || hsg[img] != hid[j]) continue;
          sigma[j] = img, sigma[img] = j;
          unres[j] = unres[img] = 0;
          mv.push_back(j), mv.push_back(img);
          left -= 2;
          ++progress;
        }
        if (progress > 0) {
          lenient = false;
          continue;
        }
        // No unique mover: columns that look unchanged stay, those whose rows meet no unresolved
        // column first (then all of them).
        for (Index j : big) {
          if (!unres[j] || hsg[j] != hid[j]) continue;
          bool clear = true;
          for (Index e = p.A.begin(j); e < p.A.end(j) && clear && !lenient; ++e) {
            const Index r = p.A.index[e];
            for (Index k = R.begin(r); k < R.end(r) && clear; ++k) clear = R.index[k] == j || !unres[R.index[k]];
          }
          if (!clear) continue;
          unres[j] = 0;
          --left;
          ++progress;
        }
        if (progress == 0) {
          if (lenient) {
            if (dbg) std::fprintf(stderr, "ORB   swap %zu: %zu columns unresolved\n", b, left);
            ok = false;
          }
          lenient = true;
        }
      }
      for (Index j : big) unres[j] = 0;
      if (ok) {
        ok = checker->check(sigma, mv);
        if (!ok && dbg) std::fprintf(stderr, "ORB   swap %zu: exact check failed (%zu moved)\n", b, mv.size());
      }
      if (ok) {
        swaps.emplace_back();
        for (Index j : mv) swaps.back().emplace_back(j, sigma[j]);
      }
      for (Index j : mv) sigma[j] = j;
    }
    if (!ok) {
      if (dbg) std::fprintf(stderr, "ORB comp of %zu cols (%zu orbits, p %zu, %zu gens): block swap not verified\n",
                            cols.size(), orbits.size(), np, comp_gens[root].size());
      cleanup();
      continue;
    }
    // Matrix rows: the orbits that are set packing rows.
    Orbitope o;
    o.cols = static_cast<Index>(blocks.size());
    o.swaps = std::move(swaps);
    std::vector<Index> entries;
    for (size_t t : small_ids) {
      entries.clear();
      for (const auto& bl : blocks) entries.push_back(bl[t]);
      if (!is_packing(p, R, entries, mark)) continue;
      o.var.insert(o.var.end(), entries.begin(), entries.end());
      ++o.rows;
    }
    if (dbg) std::fprintf(stderr, "ORB comp of %zu cols: %zu x %zu blocks, %d packing rows\n", cols.size(),
                          orbits.size(), blocks.size(), o.rows);
    cleanup();
    if (o.rows == 0) continue;
    if (used_gen) {
      for (int g : comp_gens[root]) (*used_gen)[g] = 1;
    }
    out.push_back(std::move(o));
  }
  return out;
}

void orbitope_canonicalize(const Orbitope& o, std::vector<double>& x) {
  if (static_cast<Index>(o.swaps.size()) + 1 != o.cols) return;
  auto less = [&](Index c) {  // column c lexicographically smaller than column c + 1
    for (Index i = 0; i < o.rows; ++i) {
      const double a = std::round(x[o.at(i, c)]), b = std::round(x[o.at(i, c + 1)]);
      if (a != b) return a < b;
    }
    return false;
  };
  std::vector<double> y;
  for (bool changed = true; changed;) {
    changed = false;
    for (Index c = 0; c + 1 < o.cols; ++c) {
      if (!less(c)) continue;
      y = x;
      for (const auto& [j, img] : o.swaps[c]) y[img] = x[j];
      x.swap(y);
      changed = true;
    }
  }
}

bool orbitope_propagate(const Orbitope& o, std::vector<double>& lb, std::vector<double>& ub,
                        std::vector<Index>& tightened, std::vector<char>& fixed_one) {
  const Index r = o.rows, p = o.cols;
  auto zero = [&](Index j) { return ub[j] < 0.5; };
  auto one = [&](Index j) { return lb[j] > 0.5; };
  std::vector<Index> lo(p);
  // Lower bounds on the leading rows; entries above them are zero.
  Index prev = -1;
  for (Index c = 0; c < p; ++c) {
    Index i = prev >= r ? r : prev + 1;
    while (i < r && zero(o.at(i, c))) ++i;
    lo[c] = i;
    prev = i;
    for (Index k = 0; k < i; ++k) {
      const Index j = o.at(k, c);
      if (zero(j)) continue;
      if (one(j)) return false;
      ub[j] = 0.0;
      tightened.push_back(j);
      fixed_one.push_back(0);
    }
  }
  // Upper bounds on the leading rows from the ones; a column with a one needs all earlier
  // columns to lead strictly above it.
  Index next = r;
  for (Index c = p - 1; c >= 0; --c) {
    Index u = r;
    for (Index i = 0; i < r; ++i) {
      if (one(o.at(i, c))) {
        u = i;
        break;
      }
    }
    if (next < r) u = std::min(u, next - 1);
    next = u;
    if (u >= r) continue;
    if (lo[c] > u) return false;
    if (lo[c] == u) {
      const Index j = o.at(u, c);
      if (!one(j)) {
        if (zero(j)) return false;
        lb[j] = 1.0;
        tightened.push_back(j);
        fixed_one.push_back(1);
      }
    }
  }
  return true;
}

}  // namespace milp
