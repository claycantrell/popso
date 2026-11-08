// Outer approximation of convex quadratic rows (MipProblem::quad_rows) by gradient cuts.
//
// A row lo <= g(x) <= hi kept out of the LP has one finite side and g convex towards it (g = a x + x' S x with S
// PSD for a finite hi, NSD for a finite lo). At a point x^ that violates it, the tangent plane
// g(x^) + grad g(x^) (x - x^) underestimates (overestimates) g everywhere, so it gives a valid cut that x^ violates.
// A second-order cone row (x' x <= t^2, t >= 0: the quadratic is indefinite) has its tangents from the norm form
// ||x|| - t <= 0 instead (QuadRow::cone); a violated point is on the cone's outside there too.

#include <cmath>

#include "milp/cuts/cut.hpp"

namespace milp {

int quad_oa_cuts(const MipProblem& p, const double* x, double tol, std::vector<Cut>& out) {
  if (p.quad_rows.empty()) return 0;
  std::vector<double> grad(static_cast<size_t>(p.num_cols()), 0.0);
  std::vector<std::pair<Index, double>> g;
  int added = 0;
  for (const QuadRow& q : p.quad_rows) {
    if (!(q.violation(x) > tol)) continue;  // relative, as max_violation (a cone row by its norm form)
    // The tangent of the row's convex function (QuadRow::linearize: of g for a cone row).
    double rhs = 0.0;
    if (!q.linearize(x, grad, g, rhs)) continue;
    Cut c;
    c.origin = "quad";
    c.idx.reserve(g.size());
    c.val.reserve(g.size());
    for (const auto& [j, v] : g) {
      c.idx.push_back(j);
      c.val.push_back(v);
    }
    c.rhs = rhs;
    out.push_back(std::move(c));
    ++added;
  }
  return added;
}

namespace {

class QuadOaSeparator : public Separator {
 public:
  const char* name() const override { return "quad"; }
  void separate(const SepaContext& ctx, std::vector<Cut>& out) override {
    work_ = 0;
    for (const QuadRow& q : ctx.p.quad_rows) work_ += static_cast<long long>(q.idx.size() + 2 * q.qc.size());
    quad_oa_cuts(ctx.p, ctx.x, 1e-7, out);
  }
  long long work() const override { return work_; }

 private:
  long long work_ = 0;
};

}  // namespace

std::unique_ptr<Separator> make_quad_oa_separator() { return std::make_unique<QuadOaSeparator>(); }

}  // namespace milp
