// Tests for SparseMatrix and HVector.

#include "milp/core/sparse.hpp"
#include "test_framework.hpp"

using namespace milp;

TEST(sparse_from_triplets_basic) {
  // 3 columns x 2 rows, given out of order.
  SparseMatrix m = SparseMatrix::from_triplets(3, 2, {2, 0, 1, 0}, {1, 1, 0, 0}, {5.0, 2.0, 3.0, 1.0});
  CHECK_EQ(m.num_major, 3);
  CHECK_EQ(m.num_minor, 2);
  REQUIRE(m.nnz() == 4);
  CHECK((m.start == std::vector<Index>{0, 2, 3, 4}));
  CHECK((m.index == std::vector<Index>{0, 1, 0, 1}));
  CHECK((m.value == std::vector<double>{1.0, 2.0, 3.0, 5.0}));
}

TEST(sparse_from_triplets_duplicates_and_zeros) {
  SparseMatrix m = SparseMatrix::from_triplets(2, 3, {0, 0, 0, 1, 1, 1}, {2, 0, 2, 1, 1, 0},
                                               {1.5, 0.0, 2.5, 4.0, -4.0, 7.0});
  // col 0: row 0 explicit zero dropped, row 2 summed to 4; col 1: row 1 cancels, row 0 = 7.
  REQUIRE(m.nnz() == 2);
  CHECK((m.start == std::vector<Index>{0, 1, 2}));
  CHECK_EQ(m.index[0], 2);
  CHECK_EQ(m.value[0], 4.0);
  CHECK_EQ(m.index[1], 0);
  CHECK_EQ(m.value[1], 7.0);
}

TEST(sparse_from_triplets_empty) {
  SparseMatrix m = SparseMatrix::from_triplets(4, 0, {}, {}, {});
  CHECK_EQ(m.nnz(), 0);
  CHECK_EQ(m.start.size(), 5u);
  for (Index j = 0; j < 4; ++j) CHECK_EQ(m.begin(j), m.end(j));
}

TEST(sparse_transpose) {
  // [1 0 2; 0 3 0] in CSC (3 cols, 2 rows).
  SparseMatrix a = SparseMatrix::from_triplets(3, 2, {0, 2, 1}, {0, 0, 1}, {1.0, 2.0, 3.0});
  SparseMatrix t = a.transpose();
  CHECK_EQ(t.num_major, 2);
  CHECK_EQ(t.num_minor, 3);
  CHECK((t.start == std::vector<Index>{0, 2, 3}));
  CHECK((t.index == std::vector<Index>{0, 2, 1}));
  CHECK((t.value == std::vector<double>{1.0, 2.0, 3.0}));
  SparseMatrix tt = t.transpose();
  CHECK(tt.start == a.start);
  CHECK(tt.index == a.index);
  CHECK(tt.value == a.value);
}

TEST(sparse_append) {
  SparseMatrix m;
  m.num_minor = 3;
  const Index i0[] = {0, 2};
  const double v0[] = {1.0, -1.0};
  m.append(i0, v0, 2);
  m.append(nullptr, nullptr, 0);
  CHECK_EQ(m.num_major, 2);
  CHECK((m.start == std::vector<Index>{0, 2, 2}));
  CHECK_EQ(m.nnz(), 2);
}

TEST(hvector_clear_and_rebuild) {
  HVector h;
  h.resize(10);
  h.dense[3] = 1.0;
  h.dense[7] = 1e-20;
  h.dense[9] = -2.0;
  h.rebuild_index();
  CHECK_EQ(h.count, 2);
  CHECK_EQ(h.dense[7], 0.0);
  h.clear();
  CHECK_EQ(h.count, 0);
  for (double v : h.dense) CHECK_EQ(v, 0.0);
}
