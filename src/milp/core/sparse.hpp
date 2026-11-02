#pragma once
// Compressed sparse column/row matrices and a simple sparse vector.

#include <algorithm>
#include <vector>

#include "milp/core/types.hpp"

namespace milp {

// Compressed storage. Interpreted as CSC (major = column) or CSR (major = row)
// depending on context; `transpose()` converts between the two.
struct SparseMatrix {
  Index num_major = 0;  // number of columns for CSC, rows for CSR
  Index num_minor = 0;
  std::vector<Index> start{0};  // size num_major + 1
  std::vector<Index> index;     // minor indices
  std::vector<double> value;

  Index nnz() const { return static_cast<Index>(index.size()); }
  Index begin(Index j) const { return start[j]; }
  Index end(Index j) const { return start[j + 1]; }

  // Build from (major, minor, value) triplets; duplicates are summed and
  // explicit zeros dropped. Minor indices are sorted within each major.
  static SparseMatrix from_triplets(Index num_major, Index num_minor, const std::vector<Index>& major,
                                    const std::vector<Index>& minor, const std::vector<double>& val);

  SparseMatrix transpose() const;

  // Append a major vector (e.g. a new column in CSC).
  void append(const Index* idx, const double* val, Index len);
};

// Dense-backed sparse vector: values in a dense array plus a list of
// nonzero positions. Used heavily by the LU solves (hypersparsity).
struct HVector {
  std::vector<double> dense;
  std::vector<Index> nz;  // positions that may be nonzero
  Index count = 0;        // valid entries in nz

  void resize(Index n) {
    dense.assign(n, 0.0);
    nz.assign(n, 0);
    count = 0;
  }
  void clear() {
    if (count * 4 < static_cast<Index>(dense.size())) {
      for (Index k = 0; k < count; ++k) dense[nz[k]] = 0.0;
    } else {
      std::fill(dense.begin(), dense.end(), 0.0);
    }
    count = 0;
  }
  // Recompute nz from dense, dropping tiny values.
  void rebuild_index(double drop = 1e-14) {
    // Branch-free: the density of simplex vectors is often 20-70% (unpredictable branches).
    Index c = 0;
    double* d = dense.data();
    Index* z = nz.data();
    const Index len = static_cast<Index>(dense.size());
    for (Index i = 0; i < len; ++i) {
      const bool keep = std::abs(d[i]) > drop;
      z[c] = i;
      c += keep;
      d[i] = keep ? d[i] : 0.0;
    }
    count = c;
  }
};

}  // namespace milp
