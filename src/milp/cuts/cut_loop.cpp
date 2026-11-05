#include "milp/cuts/cut_loop.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <string>
#include <unordered_set>

#include "milp/core/interrupt.hpp"

namespace milp {

namespace {

bool cuts_disabled_by_env() {
#ifdef _MSC_VER
#pragma warning(suppress : 4996)
#endif
  const char* v = std::getenv("MILP_NO_CUTS");
  return v && *v && std::strcmp(v, "0") != 0;
}

// MILP_DEBUG_SOL: file with a known optimal point, one value per column.
std::vector<double> load_debug_solution(Index n) {
#ifdef _MSC_VER
#pragma warning(suppress : 4996)
#endif
  const char* f = std::getenv("MILP_DEBUG_SOL");
  std::vector<double> x;
  if (!f) return x;
#ifdef _MSC_VER
#pragma warning(suppress : 4996)
#endif
  std::FILE* fp = std::fopen(f, "r");
  if (!fp) return x;
  double v;
  while (std::fscanf(fp, "%lf", &v) == 1) x.push_back(v);
  std::fclose(fp);
  if (static_cast<Index>(x.size()) != n) x.clear();
  return x;
}

// The known point is only meaningful if it lies within the problem's column
// bounds (sub-MIPs restrict bounds, so the global optimum is usually outside).
void drop_if_outside(const MipProblem& p, std::vector<double>& x) {
  for (Index j = 0; j < static_cast<Index>(x.size()); ++j) {
    if (x[j] < p.col_lb[j] - 1e-6 || x[j] > p.col_ub[j] + 1e-6) {
      x.clear();
      return;
    }
  }
  // Sub-MIPs may also add rows (local branching) that exclude the point.
  if (!x.empty() && p.max_violation(x) > 1e-6) x.clear();
}

}  // namespace

CutLoopStats root_cut_loop(const MipProblem& p, DualSimplex& lp, const SolverParams& params, const Logger& log,
                           const Timer& timer, double cutoff, CutData* data, bool sub_mip,
                           std::vector<Cut>* dropped, const std::function<void(const DualSimplex&)>& on_lp,
                           const std::atomic<bool>* stop, const std::function<double()>& cutoff_fn,
                           int max_rounds) {
  if (cuts_disabled_by_env()) return CutLoopStats{};
  CutLoopOptions opt;
  opt.on_lp = on_lp;
  // Reduced-cost fixing between rounds: opt-in (MILP_CUTLOOP_RCFIX=1). Small root gains only
  // (mzzv42z, ns1830653), and Gomory cuts derived from the fixed bounds are valid below the cutoff
  // only, which root restarts (cut rows become model rows) would have to respect.
  static const bool rcfix = [] { const char* v = std::getenv("MILP_CUTLOOP_RCFIX"); return v && *v && *v != '0'; }();
  if (rcfix) opt.cutoff_fn = cutoff_fn;
  opt.stop = stop;
  opt.zerohalf = !sub_mip;
  if (const char* v = std::getenv("MILP_MCF")) opt.mcf = *v != '0';
  opt.tol = params.tol;
  // Spend at most a fifth of the remaining time on cuts, not counting the initial LP solve
  // (needed anyway; on large models it can take most of that fifth).
  if (is_finite(params.time_limit)) {
    opt.time_limit = timer.elapsed() + 0.2 * (params.time_limit - timer.elapsed());
    opt.time_share = 0.2;
    opt.total_time_limit = params.time_limit;
  }
  // Flat rounds on large models with expensive separation (neos-960392: 12 rounds of ~0.5-1 s, bound
  // never moved) are capped by time; the bin packing like models they are for have cheap rounds.
  if (!params.deterministic) opt.flat_time_factor = 3.0;
  if (const char* v = std::getenv("MILP_FLAT_TIME")) opt.flat_time_factor = std::atof(v);
  // MILP_FLAT_WORK = f (default 3; 0: the time cap above): that cap in work units, in both modes (a
  // wall-clock cap moves with the machine's load: neos-3381206-awhea's 7 flat rounds before its bound moves
  // take 0.4-0.8 s unloaded and more than the 1 s floor loaded, and then the tree cannot close 415 -> 453;
  // --deterministic had no cap: neos-960392 15 rounds, now 4). MILP_FLAT_WORK_MIN: the floor.
  opt.flat_work_factor = 3.0;
  if (const char* v = std::getenv("MILP_FLAT_WORK")) opt.flat_work_factor = std::atof(v);
  if (const char* v = std::getenv("MILP_FLAT_WORK_MIN")) opt.flat_work_min = std::atof(v);
  if (const char* v = std::getenv("MILP_CHEAP_ROUNDS")) opt.max_rounds_cheap = std::atoi(v);
  if (const char* v = std::getenv("MILP_MAX_CUTS")) opt.max_cuts_per_round = std::atoi(v);
  if (const char* v = std::getenv("MILP_CHEAP_RECENT_REL")) opt.cheap_recent_rel = std::atof(v);
  // MILP_CUT_AGE = k: slack cut rows stay for k further slack re-solves (0: removed at the first).
  if (const char* v = std::getenv("MILP_CUT_AGE")) opt.cut_age = std::atoi(v);
  if (max_rounds >= 0) {
    opt.max_rounds = std::min(opt.max_rounds, max_rounds);
    opt.max_rounds_cheap = 0;
  }
  return cut_loop(p, lp, opt, log, timer, cutoff, nullptr, data, dropped);
}

CutLoopStats cut_loop(const MipProblem& p, DualSimplex& lp, const CutLoopOptions& opt_in, const Logger& log,
                      const Timer& timer, double cutoff, std::vector<Cut>* all_cuts, CutData* shared_data,
                      std::vector<Cut>* dropped_out) {
  CutLoopOptions opt = opt_in;
  CutLoopStats stats;
  const Index n = p.num_cols(), m0 = p.num_rows();
  LpOptions lpo;
  lpo.tol = opt.tol;
  lpo.objective_cutoff = cutoff;
  auto solve = [&]() {
    lpo.time_limit = std::max(0.0, opt.time_limit - timer.elapsed());
    return lp.solve(lpo);
  };
  const long long iters_start = lp.iterations();
  const double t_lp0 = timer.elapsed();
  if (opt.time_share > 0.0 && is_finite(opt.total_time_limit)) {
    // The initial solve gets the whole remaining time; the share counts from its end.
    lpo.time_limit = std::max(0.0, opt.total_time_limit - timer.elapsed());
    if (lp.solve(lpo) != LpStatus::Optimal) return stats;
    const double now = timer.elapsed();
    opt.time_limit = std::max(opt.time_limit, now + opt.time_share * (opt.total_time_limit - now));
  } else if (solve() != LpStatus::Optimal) {
    return stats;
  }
  const long long iters0 = lp.iterations();
  if (opt.on_lp) opt.on_lp(lp);
  const double t_lp_end = timer.elapsed(), initial_lp_time = t_lp_end - t_lp0;
  log.log(2, "  initial LP: %lld iterations, %.2fs (%d rows, %d cols)\n", iters0 - iters_start, initial_lp_time,
          static_cast<int>(lp.num_rows()), static_cast<int>(lp.num_cols()));
  // Work of the loop after the initial solve (flat_work_factor): the cut probing budget (its time counts for
  // the time cap too), per round the re-solve's iterations x (rows + columns) plus the separation, estimated
  // as Gomory's: a tableau row (a BTRAN and a pass over the row matrix) for each fractional integer column,
  // at most 500 (neos-960392: a flat round costs ~0.4x the initial solve, mostly Gomory).
  const double initial_work = static_cast<double>(iters0 - iters_start) * static_cast<double>(lp.num_rows() + n);
  double loop_work = 0.0;
  double lp_time = 0.0;
  int extra_rounds = 0;
  double obj = lp.objective();
  stats.bound_before = stats.bound_after = obj;

  std::vector<std::unique_ptr<Separator>> seps;
  if (opt.gomory) seps.push_back(make_gomory_separator());
  if (opt.cover) seps.push_back(make_cover_separator());
  if (opt.mir) seps.push_back(make_mir_separator());
  if (opt.mir && std::getenv("MILP_ROOT_TMIR")) seps.push_back(make_tableau_mir_separator());
  if (opt.mcf) seps.push_back(make_mcf_separator());
  if (opt.clique) seps.push_back(make_clique_separator());
  if (opt.implied_bound && !std::getenv("MILP_NO_IMPLBD")) seps.push_back(make_implied_bound_separator());
  if (opt.zerohalf && !std::getenv("MILP_NO_ZEROHALF")) seps.push_back(make_zerohalf_separator());
  if (const char* only = std::getenv("MILP_ONLY_SEPA")) {
    // Diagnostics: keep only the named separators (comma list, e.g. "mir,implbd").
    const std::string list = std::string(",") + only + ",";
    seps.erase(std::remove_if(seps.begin(), seps.end(),
                              [&](const std::unique_ptr<Separator>& s) {
                                return list.find(std::string(",") + s->name() + ",") == std::string::npos;
                              }),
               seps.end());
  }
  if (seps.empty()) return stats;
  CutData local_data;
  if (!shared_data) local_data.build(p);
  CutData& data = shared_data ? *shared_data : local_data;
  if ((opt.clique || opt.implied_bound) && !data.probed) {
    const long long nnz = static_cast<long long>(p.A.index.size());
    static const char* env_work = std::getenv("MILP_PROBE_WORK");
    const long long work = opt.probe_work >= 0 ? opt.probe_work
                           : env_work           ? std::atoll(env_work)
                                                : std::clamp<long long>(2000 * nnz, 1000000LL, 500000000LL);
    const double t0 = timer.elapsed();
    // At the time limit the probing stops (its work budget alone took up to 6.5 s on fiball at 1 thread: a
    // 10 s limit overran by 10 s on a loaded machine).
    data.probe(p, lp.x(), work, [&] { return timer.elapsed() > opt.total_time_limit; });
    loop_work += static_cast<double>(work);  // its budget (counted by the time cap too)
    log.log(2, "  cut probing: %zu implied bounds, %zu literal conflicts, %zu implied variable bounds, %.2fs\n",
            data.probe_imps.size(), data.conflict_lit.size() / 2, data.num_implied_vbs, timer.elapsed() - t0);
  }

  CutPoolOptions popt;
  popt.max_cuts = opt.max_cuts_per_round > 0 ? opt.max_cuts_per_round : std::clamp<int>(n / 10, 50, 500);
  // Density limit: a quarter of the columns, but at least twice the longest problem row.
  Index max_row = 0;
  {
    std::vector<Index> len(m0, 0);
    for (Index i : p.A.index) ++len[i];
    for (Index l : len) max_row = std::max(max_row, l);
  }
  popt.max_nnz = std::max<Index>({50, n / 4 + 20, 2 * max_row});
  if (const char* v = std::getenv("MILP_SEL_HYBRID")) popt.hybrid = *v && *v != '0';
  if (const char* v = std::getenv("MILP_SEL_PAR")) popt.max_parallelism_weak = std::atof(v);
  if (const char* v = std::getenv("MILP_SEL_PARGOOD")) popt.max_parallelism = std::atof(v);
  if (const char* v = std::getenv("MILP_SEL_GOOD")) popt.good_score = std::atof(v);
  if (const char* v = std::getenv("MILP_SEL_OBJW")) popt.obj_weight = std::atof(v);
  if (const char* v = std::getenv("MILP_SEL_INTW")) popt.int_weight = std::atof(v);
  // Cut rows: up to max(1000, 2 m0); sparse cuts up to twice that while their nonzeros
  // stay below twice those of the problem (degenerate models may need a cut per bin and pattern).
  const Index max_cut_rows = std::max<Index>(1000, 2 * m0);
  const long long max_cut_nnz = 2 * static_cast<long long>(p.A.index.size());
  std::vector<double> history{obj};
  const Index r0 = lp.num_rows();
  // Rows the loop may remove (slack or useless cuts): from the first cut row, or only its own.
  const Index first_removable = opt.keep_rows ? r0 : m0;
  std::vector<Cut> row_cuts;  // the cut in LP row r0 + k
  std::vector<Cut> dropped;   // cuts removed for slack, candidates for re-adding
  const size_t max_dropped = 20000;
  // MILP_CUT_TRACE=1 (diagnostics, verbosity 2): per round, the cuts re-added from the dropped list,
  // fresh cuts with the coefficients of a cut that was in the LP before (re-separated), the
  // selection's rejections and the removed rows by their age in rounds.
  static const bool trace = [] { const char* v = std::getenv("MILP_CUT_TRACE"); return v && *v && *v != '0'; }();
  std::vector<int> row_round;  // trace: round in which row_cuts[k] entered the LP
  std::unordered_set<uint64_t> seen_rows;
  long long tr_added = 0, tr_readded = 0, tr_resep = 0, tr_removed = 0;
  auto cut_hash = [](const Cut& c) {
    uint64_t h = 1469598103934665603ULL;
    auto mix = [&](uint64_t v) { h = (h ^ v) * 1099511628211ULL; };
    for (size_t k = 0; k < c.idx.size(); ++k) {
      uint64_t b;
      std::memcpy(&b, &c.val[k], sizeof b);
      mix(static_cast<uint64_t>(c.idx[k])), mix(b);
    }
    return h;
  };
  // Cut row aging (cut_age = k > 0): a cut row whose logical is basic with slack stays in the LP until
  // k + 1 re-solves in a row found it slack (SCIP / HiGHS keep slack rows for 10 LP solves). Removing
  // it at the first slack re-solve (k = 0) churns: on tr12-30 most removed rows were added in the same
  // round, and a fifth of the added cuts were re-adds or re-separations of removed ones. Aged rows give
  // way when a full round would not fit below max_cut_rows. At the end of the loop (keep_age < 0) every
  // cut row with a basic logical leaves, also at zero slack (the LP optimum stays): aging leaves many
  // degenerate rows behind (mzzv11 with k = 10: 1926 cut rows instead of 579), and the tree's LPs keep
  // the root's cut rows.
  std::vector<int> row_age;  // per row_cuts: consecutive slack re-solves
  auto remove_slack = [&](int keep_age, int round, int* age_hist) {
    std::vector<Index> slack_rows;
    std::vector<char> slack(lp.num_rows(), 0);
    for (Index i = first_removable; i < lp.num_rows(); ++i) {
      if (!lp.row_is_basic(i)) {
        if (i >= r0) row_age[i - r0] = 0;
        continue;
      }
      const double act = -lp.var_value(n + i), ub = -lp.var_lower(n + i);
      const bool is_slack = ub - act > 1e-6 * std::max(1.0, std::abs(ub)) || (keep_age < 0 && i >= r0);
      if (i >= r0) row_age[i - r0] = is_slack ? row_age[i - r0] + 1 : 0;
      if (is_slack && (i < r0 || row_age[i - r0] > keep_age)) slack_rows.push_back(i), slack[i] = 1;
    }
    if (keep_age > 0 &&
        lp.num_rows() - m0 - static_cast<Index>(slack_rows.size()) > max_cut_rows - popt.max_cuts) {
      for (Index i = r0; i < lp.num_rows(); ++i) {  // no room for a full round: all slack rows go
        if (!slack[i] && row_age[i - r0] > 0) slack_rows.push_back(i), slack[i] = 1;
      }
      std::sort(slack_rows.begin(), slack_rows.end());
    }
    lp.remove_rows(slack_rows);
    size_t w = 0;
    for (size_t k = 0; k < row_cuts.size(); ++k) {
      if (slack[r0 + k]) {
        if (trace) {
          const int age = round - row_round[k] + 1;  // removed rows by rounds in the LP: 1, 2, 3-5, 6+
          ++age_hist[age <= 1 ? 0 : age == 2 ? 1 : age <= 5 ? 2 : 3];
        }
        if (dropped.size() < max_dropped) dropped.push_back(std::move(row_cuts[k]));
      } else {
        if (w != k) row_cuts[w] = std::move(row_cuts[k]), row_age[w] = row_age[k];
        if (trace && w != k) row_round[w] = row_round[k];
        ++w;
      }
    }
    row_cuts.resize(w);
    row_age.resize(w);
    if (trace) row_round.resize(w);
    return slack_rows.size();
  };
  const double itol = opt.tol.integrality;
  std::vector<double> debug_sol = load_debug_solution(n);
  drop_if_outside(p, debug_sol);
  int rc_fixed_total = 0;  // reduced-cost fixings in the LP (below)
  double rc_cutoff = kInf;
  auto record = [&](const std::vector<Cut>& cuts) {
    if (all_cuts) all_cuts->insert(all_cuts->end(), cuts.begin(), cuts.end());
    for (const Cut& c : cuts) {
      if (debug_sol.empty()) break;
      // Reduced-cost fixings below the cutoff may exclude the debug solution legitimately.
      if (rc_fixed_total > 0 && p.objective(debug_sol) - p.offset >= rc_cutoff - 1e-9) break;
      const double viol = c.activity(debug_sol.data()) - c.rhs;
      if (viol > 1e-6 * std::max(1.0, std::abs(c.rhs))) {
        log.log(0, "DEBUG-SOL: %s cut violated by the debug solution (%g)\n", c.origin, viol);
      }
    }
  };

  // Reduced-cost fixing in the LP between rounds against the current cutoff (incumbents found
  // while the loop runs, e.g. by the FJ thread): nonbasic integers whose reduced cost exceeds the
  // slack are fixed in the LP for the rest of the loop (as SCIP does at the root), so later
  // rounds separate a tighter LP. Changing the opposite bound of a nonbasic column keeps the basis
  // optimal; the original bounds are restored at the end. Cuts derived from these bounds (Gomory)
  // are valid for solutions better than the cutoff, as the search's own reduced-cost fixings.
  std::vector<std::pair<Index, std::pair<double, double>>> rc_saved;  // column, original bounds
  std::vector<char> rc_changed(n, 0);
  auto rc_fix = [&]() {
    if (!opt.cutoff_fn) return;
    static const char* diag = std::getenv("MILP_CL_CUTOFF");  // diagnostics: a known cutoff (internal sense, no offset)
    const double co = diag ? std::min(std::atof(diag), opt.cutoff_fn()) : opt.cutoff_fn();
    const double z = lp.objective(), slack = co - z;
    if (!is_finite(co) || !(slack >= 0.0)) return;
    rc_cutoff = std::min(rc_cutoff, co);
    for (Index j = 0; j < n; ++j) {
      if (!p.is_int[j]) continue;
      const double lo = lp.col_lower(j), up = lp.col_upper(j);
      if (lo == up) continue;
      const VarStatus st = lp.status(j);
      const double d = lp.reduced_cost(j);
      double nlo = lo, nup = up;
      if (st == VarStatus::AtLower && d > 1e-7 && is_finite(lo)) {
        nup = std::min(up, lo + std::floor(slack / d + 1e-9));
      } else if (st == VarStatus::AtUpper && d < -1e-7 && is_finite(up)) {
        nlo = std::max(lo, up - std::floor(slack / -d + 1e-9));
      }
      if (nlo > lo + 0.5 || nup < up - 0.5) {
        if (!rc_changed[j]) rc_changed[j] = 1, rc_saved.push_back({j, {lo, up}});
        lp.set_col_bounds(j, nlo, nup);
        ++rc_fixed_total;
      }
    }
  };

  bool lp_optimal = true;     // the LP is at an optimum (the final slack cleanup needs its basis)
  int degenerate_rounds = 0;  // rounds that ended on a dual degenerate LP
  // Dual degeneracy: share of nonbasic structurals with zero reduced cost. On such LPs
  // cuts often only move the solution to an alternative optimum for several rounds
  // (bin packing like models with many identical bins) before the bound moves.
  auto is_dual_degenerate = [&]() {
    Index nonbasic = 0, degenerate = 0;
    for (Index j = 0; j < n; ++j) {
      if (lp.status(j) == VarStatus::Basic || lp.var_lower(j) == lp.var_upper(j)) continue;
      ++nonbasic;
      degenerate += std::abs(lp.reduced_cost(j)) <= 1e-9;
    }
    return nonbasic > 0 && degenerate >= 0.8 * nonbasic;
  };
  bool last_degenerate = is_dual_degenerate();
  const bool initial_degenerate = last_degenerate;
  log.log(2, "  initial LP %s dual degenerate\n", initial_degenerate ? "is" : "is not");
  std::vector<double> sep_time(seps.size(), 0.0);
  auto cheap = [&]() {
    return lp.iterations() - iters0 <= std::max<long long>(opt.cheap_work * (iters0 - iters_start), 20000);
  };
  for (int round = 1; round <= std::max(opt.max_rounds, opt.max_rounds_cheap); ++round) {
    if (timer.elapsed() >= opt.time_limit) break;
    if ((opt.stop && opt.stop->load()) || interrupted()) break;
    if (round > opt.max_rounds && !cheap()) break;
    rc_fix();
    const double* x = lp.x();
    bool frac = false;
    for (Index j = 0; j < n && !frac; ++j) {
      if (p.is_int[j]) frac = std::abs(x[j] - std::round(x[j])) > itol;
    }
    if (!frac) break;
    if (opt.flat_work_factor > 0.0) {
      Index nfrac = 0;
      for (Index j = 0; j < n; ++j) nfrac += p.is_int[j] && std::abs(x[j] - std::round(x[j])) > itol;
      loop_work += static_cast<double>(std::min<Index>(nfrac, 500)) *
                   static_cast<double>(lp.row_matrix().nnz() + lp.num_rows());
    }

    // Separate. Cuts dropped earlier for slack compete again when violated (on
    // degenerate LPs the solution otherwise cycles between alternative optima).
    SepaContext ctx{p, lp, x, m0, data, round, popt.max_nnz};
    // On dual degenerate LPs (cuts mostly move the solution between alternative optima, many cheap
    // rounds are typical) separators skip their expensive variants (c-MIR paths through integer
    // columns): neos-3381206-awhea's flat rounds are capped by time, and once its bound moves the
    // extra cuts derail it (root 100% -> 84% closed, solve 10 s -> TL; mzzv11 is not degenerate).
    ctx.expensive = !initial_degenerate &&
                    !(last_degenerate && obj - stats.bound_before <= 1e-9 * std::max(1.0, std::abs(obj)));
    CutPool pool(p, popt);
    std::vector<Cut> cand;
    int found[16] = {0}, raw[16] = {0}, from_pool = 0;
    for (size_t k = 0; k < dropped.size(); ++k) {
      const Cut& c = dropped[k];
      if (c.activity(x) - c.rhs <= popt.min_violation * std::max(1.0, std::abs(c.rhs))) continue;
      Cut copy = c;
      copy.pool_id = static_cast<int>(k);
      from_pool += pool.add(std::move(copy), x);
    }
    for (size_t s = 0; s < seps.size(); ++s) {
      cand.clear();
      const double ts = timer.elapsed();
      seps[s]->separate(ctx, cand);
      sep_time[s] += timer.elapsed() - ts;
      raw[s] = static_cast<int>(cand.size());
      for (Cut& c : cand) found[s] += pool.add(std::move(c), x);
    }
    const int pool_candidates = static_cast<int>(pool.size());
    std::vector<Cut> cuts = pool.select();
    {
      Index rows = lp.num_rows() - m0;
      long long nnz = 0;
      for (const Cut& c : row_cuts) nnz += static_cast<long long>(c.idx.size());
      size_t keep = 0;
      for (; keep < cuts.size(); ++keep, ++rows) {
        nnz += static_cast<long long>(cuts[keep].idx.size());
        if (rows >= max_cut_rows && (rows >= 2 * max_cut_rows || nnz > max_cut_nnz)) break;
      }
      cuts.resize(keep);
    }
    int tr_round_readd = 0, tr_round_resep = 0;
    if (trace) {
      for (const Cut& c : cuts) {
        const bool again = !seen_rows.insert(cut_hash(c)).second;
        tr_round_readd += c.pool_id >= 0;
        tr_round_resep += again && c.pool_id < 0;
      }
      tr_added += static_cast<long long>(cuts.size());
      tr_readded += tr_round_readd, tr_resep += tr_round_resep;
    }
    {
      // Selected pool cuts leave the pool (they become rows again).
      std::vector<char> taken(dropped.size(), 0);
      bool any = false;
      for (Cut& c : cuts) {
        if (c.pool_id >= 0) taken[c.pool_id] = 1, any = true;
        c.pool_id = -1;
      }
      if (any) {
        size_t w = 0;
        for (size_t k = 0; k < dropped.size(); ++k) {
          if (taken[k]) continue;
          if (w != k) dropped[w] = std::move(dropped[k]);
          ++w;
        }
        dropped.resize(w);
      }
    }
    if (cuts.empty()) {
      if (log.verbosity() >= 3) {
        for (size_t s = 0; s < seps.size(); ++s) log.log(3, "  %s: %d candidates, none accepted\n", seps[s]->name(), raw[s]);
      }
      break;
    }

    // Add as rows -inf <= a x <= rhs.
    std::vector<Index> start{0}, idx;
    std::vector<double> val, lo, up;
    for (const Cut& c : cuts) {
      idx.insert(idx.end(), c.idx.begin(), c.idx.end());
      val.insert(val.end(), c.val.begin(), c.val.end());
      start.push_back(static_cast<Index>(idx.size()));
      lo.push_back(-kInf);
      up.push_back(c.rhs);
    }
    const long long round_iters = lp.iterations();
    const Index first_new = lp.num_rows();
    const Basis basis_before = lp.get_basis();
    lp.add_rows(static_cast<Index>(cuts.size()), start, idx, val, lo, up);
    row_cuts.insert(row_cuts.end(), cuts.begin(), cuts.end());
    row_age.resize(row_cuts.size(), 0);
    if (trace) row_round.resize(row_cuts.size(), round);
    const double t_solve = timer.elapsed();
    const LpStatus st = solve();
    const double round_lp_time = timer.elapsed() - t_solve;
    lp_time += round_lp_time;
    // Expensive re-solves (mzzv11: 300-500 cuts per round took 11-24k iterations each, so only 6 rounds
    // fit in the time share; 250 cuts per round re-solve in ~1k): when a re-solve needs more than a
    // tenth of the initial solve's iterations (and > 5000), halve the cuts per round. Counted in
    // iterations, not time, so that deterministic runs stay reproducible.
    const long long resolve_iters = lp.iterations() - round_iters;
    loop_work += static_cast<double>(resolve_iters) * static_cast<double>(lp.num_rows() + n);
    if (resolve_iters > std::max<long long>(5000, (iters0 - iters_start) / 10) && popt.max_cuts > 50) {
      popt.max_cuts = std::max(50, popt.max_cuts / 2);
      log.log(2, "  slow LP re-solve (%lld iterations): at most %d cuts per round\n", resolve_iters, popt.max_cuts);
    }
    if (st == LpStatus::Infeasible || st == LpStatus::CutoffReached) {
      lp_optimal = false;
      // Valid cuts: the node is infeasible / cut off; keep the rows, B&B re-detects it.
      stats.cuts_added += static_cast<int>(cuts.size());
      record(cuts);
      stats.rounds = round;
      log.log(2, "  cut round %d: LP %s\n", round, to_string(st));
      break;
    }
    if (st != LpStatus::Optimal) {
      // Numerical trouble or limits: drop this round's cuts and restore the LP.
      std::vector<Index> rows;
      for (Index i = first_new; i < lp.num_rows(); ++i) rows.push_back(i);
      lp.remove_rows(rows);
      row_cuts.resize(row_cuts.size() - rows.size());
      row_age.resize(row_cuts.size());
      if (trace) row_round.resize(row_cuts.size());
      lp.set_basis(basis_before);  // the optimal basis of the previous round (removing nonbasic rows breaks it)
      lp_optimal = solve() == LpStatus::Optimal;
      log.log(2, "  cut round %d: LP %s, round discarded\n", round, to_string(st));
      break;
    }
    stats.rounds = round;
    stats.cuts_added += static_cast<int>(cuts.size());
    record(cuts);
    if (opt.on_lp) opt.on_lp(lp);
    obj = lp.objective();
    stats.bound_after = obj;

    // Drop cuts that are no longer binding (basic logical with positive slack), after cut_age further slack
    // re-solves (see remove_slack).
    int tr_age[4] = {0, 0, 0, 0};
    const size_t removed = remove_slack(opt.cut_age, round, tr_age);
    if (trace) {
      const CutPool::SelectStats& ss = pool.select_stats();
      tr_removed += static_cast<long long>(removed);
      log.log(2, "  churn round %2d: %d candidates (%d from the dropped list), %d too parallel, %d over the limit; "
              "added %zu (%d re-added, %d re-separated); removed %zu (age 1: %d, 2: %d, 3-5: %d, 6+: %d)\n", round,
              pool_candidates, from_pool, ss.parallel, ss.over_limit, cuts.size(), tr_round_readd, tr_round_resep,
              removed, tr_age[0], tr_age[1], tr_age[2], tr_age[3]);
    }

    const bool dual_degenerate = is_dual_degenerate();
    last_degenerate = dual_degenerate;
    degenerate_rounds += dual_degenerate;

    if (log.verbosity() >= 2) {
      char buf[400];
      int len = 0;
      for (size_t s = 0; s < seps.size(); ++s) {
        len += std::snprintf(buf + len, sizeof(buf) - len, " %s %d", seps[s]->name(), found[s]);
      }
      if (from_pool > 0) len += std::snprintf(buf + len, sizeof(buf) - len, " pool %d", from_pool);
      log.log(2, "  cut round %2d: bound %.10g, added %zu (found%s), %d cut rows, %lld LP iters (%.2fs), %s%.2fs\n", round,
              p.display_sense * (obj + p.offset), cuts.size(), buf, static_cast<int>(lp.num_rows() - m0),
              lp.iterations() - round_iters, round_lp_time, dual_degenerate ? "degenerate, " : "", timer.elapsed());
    }

    history.push_back(obj);
    const size_t h = history.size();
    if (h > static_cast<size_t>(opt.stall_rounds)) {
      const double old = history[h - 1 - opt.stall_rounds];
      if (obj - old < opt.stall_rel * std::max(1.0, std::abs(obj))) {
        // While the bound has not moved at all on a degenerate LP, keep going for a few
        // extra rounds within an LP work budget (a multiple of the initial solve).
        const bool flat = obj - stats.bound_before <= 1e-9 * std::max(1.0, std::abs(obj));
        const bool flat_time_ok =
            opt.flat_work_factor > 0.0
                ? loop_work <= std::max(opt.flat_work_min, opt.flat_work_factor * initial_work)
                : opt.flat_time_factor <= 0.0 ||
                      timer.elapsed() - t_lp_end <= std::max(1.0, opt.flat_time_factor * initial_lp_time);
        if (opt.flat_work_factor > 0.0 && flat && dual_degenerate) {
          log.log(3, "  flat round %d: loop work %.3g of max(%.3g, %.3g x initial %.3g)\n", round, loop_work,
                  opt.flat_work_min, opt.flat_work_factor, initial_work);
        }
        if (flat && dual_degenerate && !flat_time_ok) break;
        if (flat && dual_degenerate && extra_rounds < opt.degenerate_rounds &&
            lp.iterations() - iters0 <= 5 * std::max<long long>(iters0 - iters_start, 1000)) {
          ++extra_rounds;
        } else {
          // Cheap rounds. Mostly degenerate LPs: many weak rounds are typical (awhea), keep going
          // while the bound still moves. Otherwise only while the bound curve has not flattened:
          // the last 10 rounds gave at least cheap_recent of the total gain (gap instances).
          bool more;
          if (2 * degenerate_rounds >= round) {
            more = round < opt.cheap_min_rounds ||
                   obj - history[h > 20 ? h - 21 : 0] >= opt.cheap_rel * std::max(1.0, std::abs(obj));
          } else {
            const double recent = obj - (h > 10 ? history[h - 11] : stats.bound_before);
            more = recent >= opt.cheap_recent * (obj - stats.bound_before) && recent > 0.0;
            // Or while the last 10 rounds still moved the bound by cheap_recent_rel relative to its
            // magnitude (kasai, roll3000, beasleyC3: a long, steady tail that the tree needs).
            if (!more && opt.cheap_recent_rel > 0.0) {
              more = recent >= opt.cheap_recent_rel * std::max(1.0, std::abs(obj));
            }
          }
          if (!(more && cheap())) break;
        }
      }
    }
  }
  if (log.verbosity() >= 2 && stats.rounds > 0) {
    char buf[400];
    int len = 0;
    for (size_t s = 0; s < seps.size(); ++s) {
      len += std::snprintf(buf + len, sizeof(buf) - len, " %s %.2fs", seps[s]->name(), sep_time[s]);
    }
    log.log(2, "  separation time:%s; LP re-solves %.2fs\n", buf, lp_time);
    // Cut rows left in the LP by separator (compare with other solvers' "cutting planes" statistics).
    std::vector<std::pair<std::string, int>> by_origin;
    for (const Cut& c : row_cuts) {
      auto it = std::find_if(by_origin.begin(), by_origin.end(), [&](const auto& e) { return e.first == c.origin; });
      if (it == by_origin.end()) by_origin.emplace_back(c.origin, 1);
      else ++it->second;
    }
    std::string s;
    for (const auto& [name, cnt] : by_origin) s += " " + name + " " + std::to_string(cnt);
    log.log(2, "  cut rows by separator:%s\n", s.c_str());
    if (trace) {
      log.log(2, "  churn: %lld cuts added (%lld re-added from the dropped list, %lld re-separated), %lld removed for "
              "slack, %zu distinct cuts in the LP at some point, %zu in the dropped list\n", tr_added, tr_readded,
              tr_resep, tr_removed, seen_rows.size(), dropped.size());
    }
  }
  if (opt.cut_age > 0 && lp_optimal && stats.rounds > 0) {
    // Cut rows with basic logicals leave now: the tree's LPs keep only the cuts the root LP vertex rests on.
    int hist[4] = {0, 0, 0, 0};
    const size_t removed = remove_slack(-1, stats.rounds, hist);
    log.log(2, "  %zu cut rows with basic logicals removed at the end of the loop\n", removed);
  }
  if (!rc_saved.empty()) {
    // Restore the bounds (boxed columns left at the wrong bound are flipped by the next solve).
    for (const auto& [j, b] : rc_saved) lp.set_col_bounds(j, b.first, b.second);
    log.log(2, "  cut loop reduced-cost fixing: %d bound changes on %zu columns\n", rc_fixed_total, rc_saved.size());
  }
  // Cuts that did not move the bound only make node LPs larger: drop them all.
  if (stats.rounds > 0 && lp.num_rows() > first_removable &&
      stats.bound_after - stats.bound_before <= 1e-6 * std::max(1.0, std::abs(stats.bound_before))) {
    std::vector<Index> rows;
    for (Index i = first_removable; i < lp.num_rows(); ++i) rows.push_back(i);
    lp.remove_rows(rows);
    log.log(2, "  cuts did not improve the bound; removed\n");
  }
  // Removing basic rows keeps the basis optimal; re-solve cheaply so that the
  // LP state is consistent for the caller.
  if (stats.rounds > 0 || !rc_saved.empty()) solve();
  if (dropped_out) {
    for (Cut& c : dropped) c.pool_id = -1;
    *dropped_out = std::move(dropped);
  }
  return stats;
}

}  // namespace milp
