// Replays an LP dumped with MILP_DUMP_LP: warm solve from the dumped basis
// under several option combinations, compared with a cold solve.
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

#include "milp/lp/simplex.hpp"

using namespace milp;

static double rd(std::istream& in) {
  std::string s;
  in >> s;
  if (s == "inf") return kInf;
  if (s == "-inf") return -kInf;
  return std::stod(s);
}

int main(int argc, char** argv) {
  std::ifstream in(argv[1]);
  int n, m, nnz;
  in >> n >> m;
  std::vector<double> c(n), cl(n), cu(n), rl(m), ru(m);
  for (int j = 0; j < n; ++j) c[j] = rd(in), cl[j] = rd(in), cu[j] = rd(in);
  for (int i = 0; i < m; ++i) rl[i] = rd(in), ru[i] = rd(in);
  in >> nnz;
  std::vector<Index> I(nnz), J(nnz);
  std::vector<double> V(nnz);
  for (int k = 0; k < nnz; ++k) in >> I[k] >> J[k] >> V[k];
  Basis b;
  b.status.resize(n + m);
  for (int k = 0; k < n + m; ++k) {
    int v;
    in >> v;
    b.status[k] = static_cast<VarStatus>(v);
  }
  double cutoff = rd(in);
  long long itlim;
  in >> itlim;
  if (std::getenv("REPLAY_ROWSCALE")) {
    // Geometric row scaling by powers of two: r_i = 2^round(-log2(sqrt(max|a| * min|a|))).
    std::vector<double> mx(m, 0.0), mn(m, kInf), r(m, 1.0);
    for (int k = 0; k < nnz; ++k) mx[I[k]] = std::max(mx[I[k]], std::abs(V[k])), mn[I[k]] = std::min(mn[I[k]], std::abs(V[k]));
    for (int i = 0; i < m; ++i)
      if (mx[i] > 0) r[i] = std::ldexp(1.0, -static_cast<int>(std::lround(std::log2(std::sqrt(mx[i] * mn[i]))))), r[i] = std::getenv("REPLAY_UPONLY") ? std::max(1.0, r[i]) : r[i];
    for (int k = 0; k < nnz; ++k) V[k] *= r[I[k]];
    for (int i = 0; i < m; ++i) rl[i] *= r[i], ru[i] *= r[i];
    std::printf("row scaling applied\n");
  }
  SparseMatrix A = SparseMatrix::from_triplets(n, m, J, I, V);
  auto run = [&](const char* label, bool warm, bool use_cutoff, long long lim, bool perturb) {
    DualSimplex s;
    s.load(A, c, cl, cu, rl, ru);
    if (warm) s.set_basis(b);
    LpOptions o;
    o.perturb = perturb;
    if (use_cutoff) o.objective_cutoff = cutoff;
    o.iteration_limit = lim;
    LpStatus st = s.solve(o);
    // Dual feasibility of the returned basis w.r.t. the true costs.
    double dinf = 0.0;
    Index worst = -1;
    std::vector<double> y = s.row_duals();
    for (Index j = 0; j < n + m; ++j) {
      const VarStatus vs = s.status(j);
      if (vs == VarStatus::Basic) continue;
      double d;
      if (j < n) {
        d = c[j];
        for (Index k = A.begin(j); k < A.end(j); ++k) d -= A.value[k] * y[A.index[k]];
      } else {
        d = -y[j - n];
      }
      const double lo = s.var_lower(j), up = s.var_upper(j);
      if (lo == up) continue;
      double inf = 0.0;
      if (vs == VarStatus::AtLower) inf = -d;
      else if (vs == VarStatus::AtUpper) inf = d;
      else inf = std::abs(d);
      if (inf > dinf) dinf = inf, worst = j;
    }
    // Primal feasibility recomputed from x.
    std::vector<double> x = s.col_values(), ax(m, 0.0);
    double pinf = 0.0;
    for (Index j = 0; j < n; ++j) {
      pinf = std::max({pinf, cl[j] - x[j], x[j] - cu[j]});
      for (Index k = A.begin(j); k < A.end(j); ++k) ax[A.index[k]] += A.value[k] * x[j];
    }
    Index wrow = -1;
    for (Index i = 0; i < m; ++i) {
      const double v = std::max(rl[i] - ax[i], ax[i] - ru[i]);
      if (v > pinf) pinf = v, wrow = i;
    }
    std::printf("%-34s %-12s obj %.10g iters %lld  dual inf %.2g  primal inf %.2g (row %d)\n", label, to_string(st),
                s.objective(), s.iterations(), dinf, pinf, wrow);
  };
  std::printf("cutoff %.10g iteration limit %lld\n", cutoff, itlim);
  run("cold", false, false, -1, true);
  run("warm", true, false, -1, true);
  run("warm, no perturb", true, false, -1, false);
  run("warm + cutoff", true, true, -1, true);
  run("warm + cutoff + itlim (probe)", true, true, itlim, true);
  run("warm + itlim", true, false, itlim, true);
  run("warm + cutoff + itlim, no perturb", true, true, itlim, false);
}
