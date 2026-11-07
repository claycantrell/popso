#include "milp/solver/solver.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <limits>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>

#include "milp/core/interrupt.hpp"
#include "milp/core/log.hpp"
#include "milp/core/timer.hpp"
#include "milp/heuristics/fpr.hpp"
#include "milp/heuristics/lmip.hpp"
#include "milp/io/mps_reader.hpp"
#include "milp/io/solution.hpp"
#include "milp/presolve/presolve.hpp"
#include "milp/qp/miqp.hpp"
#include "milp/qp/qp.hpp"
#include "milp/qp/qp_presolve.hpp"
#include "milp/search/branch_and_bound.hpp"
#include "milp/search/orbitope.hpp"
#include "milp/search/sst.hpp"
#include "milp/search/symmetry.hpp"
#include "milp/search/mip_problem.hpp"
#include "milp/solver/duals.hpp"
#include "milp/solver/postsolve.hpp"
#include "milp/solver/square_terms.hpp"
#include "milp/solver/reformulate.hpp"
#include "milp/solver/start.hpp"

namespace milp {

const char* to_string(Status s) {
  switch (s) {
    case Status::Unknown: return "UNKNOWN";
    case Status::Optimal: return "OPTIMAL";
    case Status::Infeasible: return "INFEASIBLE";
    case Status::Unbounded: return "UNBOUNDED";
    case Status::InfeasibleOrUnbounded: return "INF_OR_UNBD";
    case Status::TimeLimit: return "TIME_LIMIT";
    case Status::NodeLimit: return "NODE_LIMIT";
    case Status::MemoryLimit: return "MEMORY_LIMIT";
    case Status::IterationLimit: return "ITERATION_LIMIT";
    case Status::Interrupted: return "INTERRUPTED";
    case Status::NumericalError: return "NUMERIC";
    case Status::NonConvex: return "NONCONVEX";
  }
  return "UNKNOWN";
}

const char* to_string(LpMethod m) {
  switch (m) {
    case LpMethod::Auto: return "auto";
    case LpMethod::Simplex: return "simplex";
    case LpMethod::Barrier: return "barrier";
  }
  return "auto";
}

namespace {

// Branch-and-bound, optionally as a race of differently seeded searches that
// share incumbents; the first search to finish stops the others.
MipResult search(const MipProblem& q, const SolverParams& params, const Logger& log, const Timer& timer,
                 const SearchOptions& opts) {
  const int threads = params.threads > 0 ? params.threads : static_cast<int>(std::thread::hardware_concurrency());
  const int K = std::min(params.racing, std::max(1, threads));
  if (K <= 1 || params.deterministic) return branch_and_bound(q, params, log, timer, opts);
  log.log(1, "Racing %d searches with %d threads each\n", K, std::max(1, threads / K));
  std::atomic<bool> stop{false};
  SharedIncumbent shared;
  std::atomic<int> winner{-1};
  std::vector<MipResult> res(K);
  std::vector<std::thread> pool;
  for (int k = 0; k < K; ++k) {
    pool.emplace_back(inherit_interrupt([&, k] {
      SolverParams pk = params;
      pk.threads = std::max(1, threads / K);
      pk.seed = params.seed + 7919u * static_cast<unsigned>(k);
      pk.verbosity = k == 0 ? params.verbosity : 0;
      Logger lk(pk.verbosity);
      SearchOptions ok = opts;
      ok.stop_flag = &stop;
      ok.shared = &shared;
      ok.feasibility_jump = opts.feasibility_jump && k == 0;
      res[k] = branch_and_bound(q, pk, k == 0 ? log : lk, timer, ok);
      if (res[k].status != Status::Interrupted) {
        int expected = -1;
        if (winner.compare_exchange_strong(expected, k)) stop = true;
      }
    }));
  }
  for (auto& t : pool) t.join();
  const int w = winner.load() >= 0 ? winner.load() : 0;
  MipResult out = std::move(res[w]);
  for (int k = 0; k < K; ++k) {
    if (k == w) continue;
    out.nodes += res[k].nodes;
    out.lp_iterations += res[k].lp_iterations;
    if (res[k].has_solution && (!out.has_solution || res[k].objective < out.objective)) {
      out.has_solution = true;
      out.x = res[k].x;
      out.objective = res[k].objective;
    }
    if (out.status != Status::Optimal && out.status != Status::Infeasible) out.bound = std::max(out.bound, res[k].bound);
  }
  if (out.has_solution) out.bound = std::min(out.bound, out.objective);
  return out;
}

// Symmetry-breaking rows (search/sst.hpp) for the presolved problem q: detection, then SST rows
// and fixings for the generators not handled by packing orbitopes in the search. Returns false
// (and leaves *out alone) if nothing was added. Off with --symmetry off, MILP_NO_SYMMETRY=1 or
// MILP_SST=0. MILP_SST=1 (default): leaders with conflicts in their orbit only; 2: any orbit.
// MILP_SST_ROWS=0 keeps only the fixings; MILP_SST_MIN_FIX (default 0.1): minimum share of the
// columns fixed for the reductions to be used.
bool symmetry_rows(const MipProblem& q, const SolverParams& params, const Timer& timer, const Logger& log,
                   const SearchOptions& opts, MipProblem* out) {
  const int mode = [] {
    const char* v = std::getenv("MILP_SST");
    return v ? std::atoi(v) : 1;
  }();
  if (mode <= 0 || !params.symmetry || opts.sst_done || !opts.sub_mip_heuristics || !q.implications.empty() ||
      std::getenv("MILP_NO_SYMMETRY") != nullptr) {
    return false;
  }
  if (std::find(q.is_int.begin(), q.is_int.end(), 1) == q.is_int.end()) return false;
  const double t0 = timer.elapsed();
  SymmetryOptions so;
  // Deterministic mode: the work limit only (a binding time cap would be timing dependent).
  const double rem = std::max(0.0, params.time_limit - t0);
  so.time_limit = params.deterministic ? rem : std::max(1.0, 0.03 * rem);
  const Symmetry sym = detect_symmetry(q, so);
  if (sym.empty()) return false;
  std::vector<char> used;
  if (std::getenv("MILP_NO_ORBITOPE") == nullptr) detect_packing_orbitopes(q, sym, &used);
  SstOptions sst;
  sst.mode = mode;
  sst.seed = params.seed;
  const bool rows = [] {
    const char* v = std::getenv("MILP_SST_ROWS");
    return v ? std::atoi(v) != 0 : true;
  }();
  sst.rows = rows;
  SstResult r = sst_reductions(q, sym, used, sst);
  // Only when the fixings shrink the problem substantially: a few fixings (and their rows) mostly
  // take symmetry away from orbital branching and change the root heuristics' luck (measured on
  // miplib_wide_core: neos-950242 / neos-960392 slower, decomp2 / krka / ns1208400 neutral).
  const double min_fix = [] {
    const char* v = std::getenv("MILP_SST_MIN_FIX");
    return v ? std::atof(v) : 0.1;
  }();
  bool use = !r.empty() && static_cast<double>(r.fixed_zero.size()) >= min_fix * q.num_cols();
  if (!use && !r.empty()) log.log(2, "Symmetry: %d leaders would fix %zu columns; not used\n", r.leaders, r.fixed_zero.size());
  // Otherwise SST rows on continuous / general integer orbits (SCIP's default leader types; binaries
  // stay with orbital branching and orbitopes in the search). fhnw-binpack4-4 (infeasible 2D packing,
  // identical items): rows x_l >= x_f on the positions of interchangeable items; orbital branching on
  // the binaries alone fixed 209 variables in 96k nodes and the tree did not close in 120 s; SCIP's
  // own ablations there: SST cuts only 13 s, without symmetry handling or with integer leaders only no
  // proof in 60-120 s (280-410k nodes). By default on zero-objective problems only: on 10 core models
  // they take symmetry away from orbital branching as the integer SST rows did (neos-3381206-awhea
  // 97 leaders / 5029 rows, 1.9 s -> TL; neos-2978193-inde 17 -> 50 s; neos-4738912-atrato 45 -> 93 s;
  // n5-3, neos-1171448, neos5, seymour1, neos-827175 within +-15%; 1 seed, 120 s, 4 threads).
  // MILP_SST_CONT=1: all problems, 0: off.
  const int cont_mode = [] {
    const char* v = std::getenv("MILP_SST_CONT");
    return v ? std::atoi(v) : -1;
  }();
  bool zero_obj = true;
  for (Index j = 0; j < q.num_cols() && zero_obj; ++j) zero_obj = q.cost[j] == 0.0;
  const bool cont = cont_mode > 0 || (cont_mode < 0 && zero_obj);
  if (!use && cont) {
    SstOptions s3 = sst;
    s3.mode = 3;
    s3.rows = true;
    r = sst_reductions(q, sym, used, s3);
    use = !r.rows.empty();
  }
  if (!use) return false;
  *out = apply_sst(q, r);
  log.log(1, "Symmetry: %d generators, group order ~1e%.1f; %d leaders: %zu symmetry-breaking rows, %zu columns fixed "
          "(%.2fs%s)\n", sym.num_generators(), sym.log10_order(), r.leaders, r.rows.size(), r.fixed_zero.size(),
          timer.elapsed() - t0, r.budget ? ", budget" : "");
  return true;
}

// Restart rules (read per solve, so that tests can change them).
//  - Root restarts: the first one when at least MILP_RESTART_FRAC (10%) of the integer columns, and at
//    least MILP_RESTART_MIN (10) of them, are newly fixed globally at the root (reduced-cost fixing,
//    root strong branching, conflicts); further ones at MILP_RESTART_FRAC2 (10%), at most
//    MILP_MAX_RESTARTS (4) in all. From the second on, only if the presolve after the previous restart
//    removed at least MILP_RESTART_MINRED (5%) of the integer columns (SCIP's restartminred) and the
//    restarted roots so far plus the last one's cost again stay within MILP_RESTART_WORK (1) x the first
//    root's: LP work in deterministic mode, else seconds (the budget then at least 5% of the time limit).
//    Restarted roots from the second restart on are cheaper: at most MILP_RESTART_ROUNDS (5) cut rounds
//    (the earlier roots' cut rows are carried as rows), and root dives / RENS / root LNS only when the
//    search before the restart improved the incumbent. (The first restart at 5%, as in the 2026-11-03
//    design, cost time on quick roots: binkar10_1, neos-911970 and neos17, 8-9% fixed, restarted only
//    at 5% and were 0.4-1.4 s slower in both seeds of the core A/B.)
//  - In-tree restarts: at most MILP_MAX_TREE_RESTARTS (1) (rule in branch_and_bound.cpp).
//  - Opt-in MILP_RESTART_PC=1 / 2: restarted searches start from the pseudocosts of the search before
//    (1: as history, counts capped; 2: as estimates only; see kRestartPcDefault). Off by default: on
//    restarting core instances (deterministic, 2 threads) both modes needed ~15% more LP iterations.
// The rule before repeated restarts: MILP_MAX_RESTARTS=1 MILP_RESTART_MIN=20.
struct RestartRule {
  int max_restarts = 4, max_tree_restarts = 1;
  double frac = 0.10, frac2 = 0.10, min_removed = 0.05, work = 1.0;
  Index min_cols = 10;
  int cheap_rounds = 5;
  bool pc = true;
  static RestartRule from_env() {
    auto env = [](const char* name, double def) {
      const char* v = std::getenv(name);
      return v && *v ? std::atof(v) : def;
    };
    RestartRule r;
    r.max_restarts = static_cast<int>(env("MILP_MAX_RESTARTS", r.max_restarts));
    r.max_tree_restarts = static_cast<int>(env("MILP_MAX_TREE_RESTARTS", r.max_tree_restarts));
    r.frac = env("MILP_RESTART_FRAC", r.frac);
    r.frac2 = env("MILP_RESTART_FRAC2", r.frac2);
    r.min_cols = static_cast<Index>(env("MILP_RESTART_MIN", static_cast<double>(r.min_cols)));
    r.min_removed = env("MILP_RESTART_MINRED", r.min_removed);
    r.work = env("MILP_RESTART_WORK", r.work);
    r.cheap_rounds = static_cast<int>(env("MILP_RESTART_ROUNDS", r.cheap_rounds));
    r.pc = env("MILP_RESTART_PC", kRestartPcDefault) != 0.0;
    return r;
  }
};

Index count_ints(const MipProblem& q) {
  Index n = 0;
  for (char c : q.is_int) n += c != 0;
  return n;
}

// Solves the problem with presolve; returns false if the caller should solve
// the original problem instead (presolve gave up or postsolve failed).
bool solve_presolved(const MipProblem& p, const SolverParams& params, const Logger& log, const Timer& timer,
                     MipResult* out, const SearchOptions& opts = {}) {
  const RestartRule rule = RestartRule::from_env();
  Presolver pre;
  const PresolveStatus ps = pre.run(p, PresolveOptions::from_env(), log);
  const PresolveStats& st = pre.stats();
  log.log(1, "Presolve: %d rows, %d cols, %d nnz -> %d rows, %d cols, %d nnz (%.2fs)\n", st.rows_before,
          st.cols_before, st.nnz_before, st.rows_after, st.cols_after, st.nnz_after, st.time);
  // This presolve, for MipResult::presolve: the outermost call sets it last (restarts and the
  // symmetry-breaking rows presolve again in nested calls), so it is the first presolve's.
  PresolveSummary summary;
  summary.ran = true;
  summary.rows_before = st.rows_before, summary.cols_before = st.cols_before, summary.nnz_before = st.nnz_before;
  summary.rows = st.rows_after, summary.cols = st.cols_after, summary.nnz = st.nnz_after;
  summary.time = st.time;
  if (ps == PresolveStatus::Infeasible) {
    log.log(1, "Presolve: problem is infeasible\n");
    out->status = Status::Infeasible;
    out->bound = kInf;
    out->presolve = summary;
    return true;
  }
  // Unbounded-or-infeasible: let the full search classify it precisely.
  if (ps == PresolveStatus::Unbounded) {
    log.log(1, "Presolve: problem is unbounded or infeasible\n");
    out->presolve = summary;
    return false;
  }
  const MipProblem& q = pre.reduced();
  if (opts.sub_mip_heuristics && opts.restarts_done == 0 && opts.tree_restarts_done == 0 && !opts.sst_done &&
      local_search_bench(q)) {
    out->status = Status::Unknown;  // MILP_LS_BENCH: heuristic evaluation only
    return true;
  }
  if (const char* fb = std::getenv("MILP_FPR_BENCH")) {
    // Standalone heuristic comparison on the presolved problem (FPR vs FJ vs local MIP, same work
    // budget); no search.
    heuristic_bench(q, std::atoll(fb) > 0 ? std::atoll(fb) : 1000000000LL, params.seed, log);
    out->status = Status::Unknown;
    return true;
  }
  MipResult mr;
  std::function<void(const std::vector<double>&, double)> q_hook;  // opts.on_incumbent in q's columns
  std::vector<double> q_start;  // opts.start in q's columns
  double q_start_obj = kInf;    // its objective before the crush (SearchOptions::start_objective)
  if (q.num_cols() == 0) {
    // Everything was removed; remaining rows (if any) have no columns.
    for (Index i = 0; i < q.num_rows(); ++i) {
      if (q.row_lb[i] > params.tol.primal_feas || q.row_ub[i] < -params.tol.primal_feas) {
        out->status = Status::Infeasible;
        out->bound = kInf;
        return true;
      }
    }
    mr.status = Status::Optimal;
    mr.has_solution = true;
    mr.objective = mr.bound = mr.root_bound = q.offset;
  } else {
    SearchOptions o = opts;
    if (opts.on_incumbent) {
      // Incumbents of q are reported in p's columns.
      o.on_incumbent = [&pre, outer = opts.on_incumbent](const std::vector<double>& x, double bound) {
        outer(pre.postsolve(x), bound);
      };
    }
    q_hook = o.on_incumbent;
    if (!o.start.empty() && static_cast<Index>(o.start.size()) != p.num_cols()) {
      // A start in other columns (a caller's bug; crush would read it as this problem's point) is dropped.
      log.log(1, "WARNING: start point with %zu values for %d columns ignored\n", o.start.size(), p.num_cols());
      o.start.clear();
    }
    if (!o.start.empty()) {
      o.start_objective = p.objective(o.start);  // the start's own objective (its crushed point may differ)
      o.start = pre.crush(o.start);
    }
    q_start = o.start;
    q_start_obj = o.start_objective;
    if (!o.pc_init.empty()) o.pc_init = o.pc_init.mapped(pre.col_map());
    o.allow_restart = opts.allow_restart && opts.restarts_done < rule.max_restarts;
    o.restart_frac = opts.restarts_done == 0 ? rule.frac : rule.frac2;
    o.restart_min_cols = rule.min_cols;
    o.restart_min_removed = rule.min_removed;
    {
      // After a restart: the share of the integer columns this presolve removed (incl. the fixings).
      const Index before = count_ints(p), after = count_ints(q);
      o.last_restart_removed = before > 0 ? static_cast<double>(before - after) / static_cast<double>(before) : 1.0;
    }
    if (!o.debug_point.empty()) {
      o.debug_point = pre.crush(o.debug_point);
      if (q.max_violation(o.debug_point) > 1e-6) {
        log.log(0, "DEBUG-SOL: the debug point is cut off by presolve (violation %.3g)\n", q.max_violation(o.debug_point));
        o.debug_point.clear();
      }
    }
    MipProblem qs;
    if (symmetry_rows(q, params, timer, log, opts, &qs)) {
      // Solved as a presolved problem of its own (presolve uses the fixings and new rows); its
      // restarts happen inside, and its columns are q's.
      o.sst_done = true;
      if (!o.debug_point.empty() && qs.max_violation(o.debug_point) > 1e-6) {
        log.log(0, "DEBUG-SOL: the debug point is cut off by the symmetry-breaking rows (violation %.3g)\n",
                qs.max_violation(o.debug_point));
        o.debug_point.clear();
      }
      if (!solve_presolved(qs, params, log, timer, &mr, o)) mr = search(qs, params, log, timer, o);
    } else {
      o.sst_done = true;
      mr = search(q, params, log, timer, o);
    }
  }
  if (mr.restart) {
    // Root restart: tightened bounds + root cuts, presolved again; the
    // incumbent found so far is the cutoff and the fallback answer.
    const double t_restart = timer.elapsed();
    MipProblem q2 = q;
    q2.col_lb = mr.restart_lb;
    q2.col_ub = mr.restart_ub;
    if (mr.has_solution) {
      // Keep the incumbent feasible (the fixings only matter for better points) so the
      // restarted search can start from it (RINS, guided dives, local branching).
      for (Index j = 0; j < q.num_cols(); ++j) {
        q2.col_lb[j] = std::min(q2.col_lb[j], mr.x[j]);
        q2.col_ub[j] = std::max(q2.col_ub[j], mr.x[j]);
      }
    }
    const Index m0 = q.num_rows(), extra = static_cast<Index>(mr.cut_lb.size());
    if (extra > 0) {
      std::vector<Index> J, I;
      std::vector<double> V;
      for (Index j = 0; j < q.num_cols(); ++j)
        for (Index k = q.A.begin(j); k < q.A.end(j); ++k) J.push_back(j), I.push_back(q.A.index[k]), V.push_back(q.A.value[k]);
      for (Index r = 0; r < extra; ++r)
        for (Index k = mr.cut_start[r]; k < mr.cut_start[r + 1]; ++k)
          J.push_back(mr.cut_index[k]), I.push_back(m0 + r), V.push_back(mr.cut_value[k]);
      q2.A = SparseMatrix::from_triplets(q.num_cols(), m0 + extra, J, I, V);
      q2.row_lb.insert(q2.row_lb.end(), mr.cut_lb.begin(), mr.cut_lb.end());
      q2.row_ub.insert(q2.row_ub.end(), mr.cut_ub.begin(), mr.cut_ub.end());
    }
    SearchOptions o2 = opts;
    o2.on_incumbent = q_hook;  // q2 has q's columns
    o2.sst_done = true;  // symmetry-breaking rows only before the first search
    // Root and in-tree restarts are counted separately (MILP_MAX_RESTARTS / MILP_MAX_TREE_RESTARTS).
    // The estimation-driven tree restart: once per solve, apart from the fixing-driven ones.
    o2.restarts_done = opts.restarts_done + (mr.tree_restart ? 0 : 1);
    o2.tree_restarts_done = opts.tree_restarts_done + (mr.tree_restart ? 1 : 0);
    o2.est_restarts_done = opts.est_restarts_done + (mr.est_restart ? 1 : 0);
    o2.allow_restart = o2.restarts_done < rule.max_restarts;
    o2.allow_tree_restart =
        opts.allow_tree_restart && o2.tree_restarts_done - o2.est_restarts_done < rule.max_tree_restarts;
    o2.allow_est_restart = opts.allow_est_restart && !mr.est_restart;
    o2.restart_frac = o2.restarts_done == 0 ? rule.frac : rule.frac2;
    o2.restart_min_cols = rule.min_cols;
    o2.restart_min_removed = rule.min_removed;
    {
      // Budget of the restarted roots (see RestartRule): the first root's cost, less the restarted roots.
      const bool first_root = opts.restarts_done == 0 && opts.tree_restarts_done == 0;
      if (first_root) {
        o2.restart_work_left = rule.work * mr.root_lp_work;
        o2.restart_time_left =
            rule.work * (is_finite(params.time_limit) ? std::max(mr.root_time, 0.05 * params.time_limit) : mr.root_time);
      } else {
        o2.restart_work_left = opts.restart_work_left - mr.root_lp_work;
        o2.restart_time_left = opts.restart_time_left - mr.root_time;
      }
    }
    o2.restart_t0 = t_restart;
    o2.root_cut_rounds = o2.restarts_done >= 2 ? rule.cheap_rounds : -1;
    o2.root_heuristics = o2.restarts_done < 2 || mr.improved;
    o2.pc_init = rule.pc ? std::move(mr.pc) : PseudocostSnapshot{};
    // The restarted problem has q's columns: the search's incumbent, else the start it was given (o2 is a copy of
    // opts, whose start is in p's columns: drayage-25-23 on x86 restarted a second time from a search that had
    // found nothing better than its start, the stale start was crushed through the wrong presolve, the
    // restarted root was pruned against it, and the incumbent 101344.27 came back as optimal; the optimum is
    // 101282.65).
    o2.start = mr.has_solution ? mr.x : q_start;
    o2.start_objective = mr.has_solution ? q.objective(mr.x) : q_start_obj;
    o2.debug_point = mr.debug_point;
    o2.feasibility_jump = false;  // it already ran
    // After a restart from a quick root (set packing/covering: a few hundred LP iterations on a few hundred rows), a second
    // zero-half loop mostly costs time (setpack/setcover medium 0.1 -> 0.2 s); MILP_RESTART_ZEROHALF=1 keeps it.
    static const bool restart_zh = std::getenv("MILP_RESTART_ZEROHALF") != nullptr;
    o2.no_zerohalf = !restart_zh && mr.root_lp_work < 1e6;
    o2.root_lns = !mr.root_lns_ran || mr.root_lns_improved;
    MipResult mr2;
    if (!solve_presolved(q2, params, log, timer, &mr2, o2)) mr2 = search(q2, params, log, timer, o2);
    mr2.nodes += mr.nodes;
    mr2.lp_iterations += mr.lp_iterations;
    mr2.root_bound = mr.root_bound;
    mr2.restarts += mr.restarts + (mr.tree_restart ? 0 : 1);
    mr2.tree_restarts += mr.tree_restarts + (mr.tree_restart ? 1 : 0);
    if (!mr2.has_solution && mr.has_solution) {
      // Nothing better than the incumbent exists (or was found in time).
      mr2.has_solution = true;
      mr2.x = std::move(mr.x);
      mr2.objective = mr.objective;
      if (mr2.status == Status::Infeasible) {
        mr2.status = Status::Optimal;
        mr2.bound = mr.objective;
      }
    }
    // A restarted search stopped at a limit (time, interrupt) may not have re-solved its root yet;
    // the bound before the restart is still valid (the restart only adds valid reductions).
    if (mr2.status != Status::Optimal && mr2.status != Status::Infeasible && is_finite(mr.bound))
      mr2.bound = std::max(mr2.bound, mr.bound);
    mr2.bound = std::min(mr2.bound, mr2.objective);
    mr = std::move(mr2);
  }
  if (!postsolve_result(p, pre, params, timer, log, &mr)) {
    out->bound = mr.bound;  // still a valid bound for the caller's re-solve
    out->presolve = summary;
    return false;
  }
  *out = std::move(mr);
  out->presolve = summary;
  return true;
}

// Reduced problem as a user model (columns c<j>, rows r<i>; minimization).
// With orig_names / col_map, columns keep their original names (diagnostics: diffs against other presolvers).
Model reduced_model(const MipProblem& q, const std::vector<std::string>* orig_names = nullptr,
                    const std::vector<Index>* col_map = nullptr) {
  Model m;
  m.name = "presolved";
  m.obj = q.cost;
  m.obj_offset = q.offset;
  m.col_lb = q.col_lb;
  m.col_ub = q.col_ub;
  m.row_lb = q.row_lb;
  m.row_ub = q.row_ub;
  m.A = q.A;
  for (Index j = 0; j < q.num_cols(); ++j) {
    m.col_type.push_back(q.is_int[j] ? VarType::Integer : VarType::Continuous);
    if (orig_names && col_map && (*col_map)[j] < static_cast<Index>(orig_names->size())) {
      m.col_names.push_back((*orig_names)[(*col_map)[j]]);
    } else {
      m.col_names.push_back("c" + std::to_string(j));
    }
  }
  for (Index i = 0; i < q.num_rows(); ++i) m.row_names.push_back("r" + std::to_string(i));
  return m;
}

// Presolve tooling (no search): MILP_PRESOLVE_ONLY=1 just presolves; MILP_PRESOLVE_DUMP=<file.mps> writes the
// reduced problem; MILP_PRESOLVE_CHECK=<file.sol> postsolves a solution of that dump and reports its violation
// and objective; MILP_PRESOLVE_ONLY=symmetry also reports symmetry detection. Returns false if no such variable is set.
bool presolve_tool(const Model& model, const MipProblem& p, const Logger& log, SolveResult* r) {
  const char* only = std::getenv("MILP_PRESOLVE_ONLY");
  const char* dump = std::getenv("MILP_PRESOLVE_DUMP");
  const char* check = std::getenv("MILP_PRESOLVE_CHECK");
  if (!only && !dump && !check) return false;
  Presolver pre;
  const PresolveStatus ps = pre.run(p, PresolveOptions::from_env(), log);
  const PresolveStats& st = pre.stats();
  log.log(1, "Presolve: %d rows, %d cols, %d nnz -> %d rows, %d cols, %d nnz (%.2fs)\n", st.rows_before,
          st.cols_before, st.nnz_before, st.rows_after, st.cols_after, st.nnz_after, st.time);
  r->status = ps == PresolveStatus::Infeasible ? Status::Infeasible
              : ps == PresolveStatus::Unbounded ? Status::InfeasibleOrUnbounded : Status::Unknown;
  if (ps != PresolveStatus::Reduced) return true;
  if (only && std::string(only) == "symmetry") {
    const Symmetry sym = detect_symmetry(pre.reduced());
    long long support = 0;
    for (const auto& g : sym.gens) support += static_cast<long long>(g.size());
    log.log(1, "Symmetry: %d generators, %d columns moved, support %lld, order ~1e%.1f, work %lld, %.3fs%s\n",
            sym.num_generators(), sym.moved_cols, support, sym.log10_order(), sym.work, sym.seconds,
            sym.aborted ? " (aborted)" : "");
  }
  const bool keep_names = std::getenv("MILP_PRESOLVE_NAMES") != nullptr;  // original column names in the dump
  const Model red = keep_names ? reduced_model(pre.reduced(), &model.col_names, &pre.col_map()) : reduced_model(pre.reduced());
  if (dump) write_mps(red, dump);
  if (check) {
    std::vector<double> xr = read_solution(check, red);
    for (double& v : xr) {
      if (v != v) v = 0.0;
    }
    const std::vector<double> x = pre.postsolve(xr);
    const double viol = p.max_violation(x);
    std::printf("Postsolve check: violation %.3g, objective %.17g (reduced objective %.17g)\n", viol,
                model.objective_value(x), red.objective_value(xr));
    r->has_solution = true;
    r->x = x;
    r->objective = model.objective_value(x);
  }
  return true;
}

// Reports new incumbents to SolverParams::on_incumbent: in the user model's columns (the first
// num_cols of the solved problem's), serialized, only strictly improving ones. A callback that asks
// to stop (or throws) sets the interrupt flag.
class IncumbentReporter {
 public:
  IncumbentReporter(const Model& model, const SolverParams& params)
      : model_(model), sense_(model.sense == ObjSense::Maximize ? -1.0 : 1.0), fn_(params.on_incumbent),
        flag_(params.interrupt) {}
  void report(const std::vector<double>& x_in, double bound) {
    std::vector<double> x(x_in);
    x.resize(static_cast<size_t>(model_.num_cols()));
    const double obj = model_.objective_value(x);
    std::lock_guard<std::mutex> lk(mu_);
    if (stopped_ || !(sense_ * obj < best_)) return;
    best_ = sense_ * obj;
    bool stop = true;
    try {
      stop = fn_(x, obj, sense_ * std::min(bound, sense_ * obj));
    } catch (...) {
    }
    if (stop) {
      stopped_ = true;  // no calls after a stop request
      if (flag_) flag_->store(true);
    }
  }

 private:
  const Model& model_;
  const double sense_;
  const std::function<bool(const std::vector<double>&, double, double)>& fn_;
  std::atomic<bool>* flag_;
  std::mutex mu_;
  double best_ = kInf;
  bool stopped_ = false;
};

// Collects SolveResult::incumbent_trace (SolverParams::trace_incumbents) from the searches'
// SearchOptions::on_trace / MiqpOptions::on_trace: only incumbents better than every one before by
// more than a relative 1e-9 (a restarted search accepts its start again, possibly with an objective
// that differs in the last bits after the new presolve; racing searches report their own).
class IncumbentTracer {
 public:
  explicit IncumbentTracer(ObjSense sense) : sense_(sense == ObjSense::Maximize ? -1.0 : 1.0) {}
  void record(double objective, const char* source, double time) {
    const double v = sense_ * objective;
    std::lock_guard<std::mutex> lk(mu_);
    const double tol = is_finite(best_) ? 1e-9 * std::max(1.0, std::abs(best_)) : 0.0;
    if (!(v < best_ - tol)) return;
    best_ = v;
    const char* src = source ? source : "";
    if (!events_.empty()) {
      IncumbentEvent& e = events_.back();
      if (e.source == src && time - e.last_time <= std::max(0.01, 0.01 * e.last_time)) {
        ++e.count;
        e.last_time = time;
        e.last_objective = objective;
        return;
      }
    }
    IncumbentEvent e;
    e.time = e.last_time = time;
    e.objective = e.last_objective = objective;
    e.source = src;
    events_.push_back(std::move(e));
  }
  std::function<void(double, const char*, double)> hook() {
    return [this](double objective, const char* source, double time) { record(objective, source, time); };
  }
  std::vector<IncumbentEvent> take() {
    std::lock_guard<std::mutex> lk(mu_);
    return std::move(events_);
  }

 private:
  const double sense_;
  std::mutex mu_;
  double best_ = kInf;
  std::vector<IncumbentEvent> events_;
};

// LP duals of a pure LP (the solution is optimal) or the fixed-MIP duals of a MILP
// (SolverParams::duals >= 2), in the model's sense, with the basis they come from (solver/duals.hpp).
void add_duals(const MipProblem& p, const SolverParams& params, const Timer& timer, const Logger& log, SolveResult& r) {
  const bool mip = std::find(p.is_int.begin(), p.is_int.end(), 1) != p.is_int.end();
  if (params.duals < (mip ? 2 : 1)) return;
  if (!mip && r.status != Status::Optimal) return;
  std::vector<double> lb, ub, x = r.x;
  if (mip) {
    lb = p.col_lb;
    ub = p.col_ub;
    for (Index j = 0; j < p.num_cols(); ++j) {
      if (!p.is_int[j]) continue;
      x[j] = std::round(x[j]);
      lb[j] = ub[j] = std::clamp(x[j], p.col_lb[j], p.col_ub[j]);
    }
  }
  const double t0 = timer.elapsed();
  // Post-processing: at least a second even when the time limit is (nearly) used up.
  const LpDualResult d = recover_lp_duals(p, lb, ub, x, std::max(1.0, params.time_limit - t0));
  log.log(1, "Duals%s: %s, %lld simplex iterations (%.2fs)\n", mip ? " (integers fixed)" : "",
          d.ok ? "recovered" : "failed", d.iterations, timer.elapsed() - t0);
  if (!d.ok) return;
  // The basis the duals come from must reproduce the solution's objective (an LP: the same optimum;
  // fixed-MIP: the incumbent is feasible for the fixed LP, so its optimum is no worse). Otherwise
  // the recovery went wrong (numerics, a point that is not optimal) and the duals would be stale.
  const double obj = p.objective(x);
  const DualObjectiveCheck check = check_dual_objective(d.objective, obj, mip);
  if (check == DualObjectiveCheck::Mismatch) {
    log.log(1, "Duals: dropped (their basis has objective %.10g, the solution %.10g)\n", p.display_sense * d.objective,
            p.display_sense * obj);
    return;
  }
  if (check == DualObjectiveCheck::BetterThanSolution) {
    log.log(1, "Duals: the fixed LP's optimum %.10g is better than the incumbent's %.10g (its continuous part is not "
            "optimal): the duals are those of the fixed LP's optimum\n", p.display_sense * d.objective,
            p.display_sense * obj);
  }
  const double sense = p.display_sense;
  r.row_dual = d.row_dual;
  r.reduced_cost = d.reduced_cost;
  for (double& v : r.row_dual) v *= sense;
  for (double& v : r.reduced_cost) v *= sense;
  r.col_basis.resize(d.col_status.size());
  r.row_basis.resize(d.row_status.size());
  for (size_t j = 0; j < d.col_status.size(); ++j) r.col_basis[j] = static_cast<int8_t>(d.col_status[j]);
  for (size_t i = 0; i < d.row_status.size(); ++i) r.row_basis[i] = static_cast<int8_t>(d.row_status[i]);
  r.has_duals = true;
  r.fixed_mip_duals = mip;
  r.stats.time = timer.elapsed();
}

// Solves a plain MILP; `implications` (logical constraints without a big-M form) go into the
// searched problem, and presolve is skipped when there are any.
SolveResult solve_plain(const Model& model, const std::vector<Implication>& implications, const SolverParams& params_in,
                        const std::vector<double>& start, const Timer& timer, const Logger& log,
                        IncumbentReporter* reporter, IncumbentTracer* tracer, bool allow_duals) {
  SolveResult r;
  MipProblem p;
  if (!MipProblem::from_model(model, &p)) {
    r.status = Status::Infeasible;
    return r;
  }
  p.implications = implications;
  SolverParams params = params_in;
  if (!implications.empty() && params.presolve) {
    log.log(1, "Presolve: skipped (%zu logical implications)\n", implications.size());
    params.presolve = false;
  }
  const double sense = model.sense == ObjSense::Maximize ? -1.0 : 1.0;
  if (params.presolve && presolve_tool(model, p, log, &r)) return r;
  MipResult mr;
  SearchOptions opts;
  opts.allow_restart = true;
  opts.allow_tree_restart = true;
  opts.allow_est_restart = true;
  if (reporter) {
    opts.on_incumbent = [reporter](const std::vector<double>& x, double bound) { reporter->report(x, bound); };
  }
  if (tracer) opts.on_trace = tracer->hook();
  if (const char* v = std::getenv("MILP_CUTOFF")) opts.cutoff = sense * std::atof(v);  // diagnostics: known bound
  if (const char* f = std::getenv("MILP_DEBUG_SOL_ORIG")) {
    // Debugging: a known optimal point of the user model (one value per column), followed through
    // presolve and restarts (MILP_DEBUG_SOL is in the columns of the searched problem instead).
    if (std::FILE* fp = std::fopen(f, "r")) {
      double v;
      while (std::fscanf(fp, "%lf", &v) == 1) opts.debug_point.push_back(v);
      std::fclose(fp);
    }
    if (static_cast<Index>(opts.debug_point.size()) != p.num_cols()) opts.debug_point.clear();
  }
  // User start: made feasible here (original columns), then followed through presolve as the first incumbent.
  std::vector<double> start_x;
  if (!start.empty()) {
    StartOutcome so = prepare_start(p, start, params, timer);
    log.log(1, "%s\n", so.message.c_str());
    start_x = std::move(so.x);
    opts.start = start_x;
  }
  if (!params.presolve || !solve_presolved(p, params, log, timer, &mr, opts)) {
    // A failed presolved solve may still leave a valid bound (postsolve failure).
    const double presolved_bound = params.presolve && is_finite(mr.bound) ? mr.bound : -kInf;
    const PresolveSummary presolved = mr.presolve;  // a presolve that gave up still ran
    SearchOptions plain;
    plain.start = start_x;
    plain.on_incumbent = opts.on_incumbent;
    plain.on_trace = opts.on_trace;
    mr = branch_and_bound(p, params, log, timer, plain);
    mr.presolve = presolved;
    if (presolved_bound > mr.bound && !(mr.has_solution && presolved_bound > mr.objective)) {
      mr.bound = presolved_bound;
    }
  }
  if (!start_x.empty()) {
    // The start is feasible; the search may have lost it (presolve reductions can cut it off, and its
    // objective then only served as a cutoff). Nothing better found means it is the answer.
    const double so = p.objective(start_x);
    if (!mr.has_solution || so < mr.objective) {
      mr.has_solution = true;
      mr.objective = so;
      mr.x = start_x;
      if (mr.status == Status::Infeasible) {
        mr.status = Status::Optimal;
        mr.bound = so;
      } else if (mr.status == Status::InfeasibleOrUnbounded) {
        mr.status = Status::Unbounded;
      }
    }
    if (mr.status != Status::Unbounded) mr.bound = std::min(mr.bound, mr.objective);
  }
  r.status = mr.status;
  r.has_solution = mr.has_solution;
  if (mr.has_solution) {
    r.x = std::move(mr.x);
    r.objective = model.objective_value(r.x);
  }
  // Internal bounds are for minimization; flip back for maximization.
  r.bound = sense * mr.bound;
  r.stats.time = timer.elapsed();
  r.stats.nodes = mr.nodes;
  r.stats.lp_iterations = mr.lp_iterations;
  r.stats.root_bound = sense * mr.root_bound;
  r.stats.restarts = mr.restarts;
  r.stats.tree_restarts = mr.tree_restarts;
  r.stats.presolve = mr.presolve;
  if (allow_duals && implications.empty() && r.has_solution && !interrupted()) add_duals(p, params, timer, log, r);
  return r;
}

// Values within the tolerance of zero that the logical constraints require to be zero (a semi-continuous
// column below its lower bound, SOS members outside the nonzero pattern) are set to exactly 0, if that
// does not make the point less feasible.
void snap_logical(const Model& model, std::vector<double>& x) {
  constexpr double kSnap = 1e-5;
  std::vector<double> y = x;
  for (Index j = 0; j < model.num_cols(); ++j) {
    if (!is_semi_type(model.col_type[j]) || std::abs(y[j]) > kSnap) continue;
    if (y[j] < model.col_lb[j] || y[j] > model.col_ub[j]) y[j] = 0.0;
  }
  for (const SosConstraint& s : model.sos) {
    // Keep the largest member (SOS1) or the pair of neighbours with the largest values (SOS2).
    const size_t k = s.cols.size();
    std::vector<size_t> ord(k);
    for (size_t t = 0; t < k; ++t) ord[t] = t;
    std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return s.weights[a] < s.weights[b]; });
    size_t keep = 0;
    double best = -1.0;
    const size_t span = s.type == 2 ? 2 : 1;
    for (size_t t = 0; t + span <= k; ++t) {
      double v = 0.0;
      for (size_t u = t; u < t + span; ++u) v = std::max(v, std::abs(y[s.cols[ord[u]]]));
      if (v > best) best = v, keep = t;
    }
    for (size_t t = 0; t < k; ++t) {
      const Index j = s.cols[ord[t]];
      if ((t < keep || t >= keep + span) && std::abs(y[j]) <= kSnap) y[j] = 0.0;
    }
  }
  if (y != x && model.violation(y).max() <= std::max(1e-9, model.violation(x).max())) x = std::move(y);
}

SolveResult solve_quadratic(const Model& model, const SolverParams& params, const std::vector<double>& start,
                            const Timer& timer, const Logger& log, IncumbentTracer* tracer);

// A separable convex objective on integer columns (square_terms.hpp) as a MILP: tangent rows of each square
// around the continuous relaxation's point, the MILP solved by the MILP pipeline, its solution scored with
// the exact objective, and tangent windows added where it falls outside them, until the solution is exact.
// Each MILP is a relaxation of the model, so its bound is the model's. MILP_SQ_WINDOW: the window (25);
// MILP_SQ_ROUNDS: at most that many MILPs (10); MILP_SQ_CENTER=0: no relaxation (centers at 0). A user's integer QP (141637 rows, 265 squared general
// integers): the MIQP search had no incumbent after 180 s (its root dive solves a barrier QP per step); as
// a MILP, OPTIMAL at the root in 4.6 s (root LP 9: the tangents give the convex hull of x^2 on the integers).
SolveResult solve_square_milp(const Model& model, SquareTerms& sq, const SolverParams& params,
                              const std::vector<double>& start, const Timer& timer, const Logger& log,
                              IncumbentTracer* tracer) {
  const double sense = model.sense == ObjSense::Maximize ? -1.0 : 1.0;
  const Index n = model.num_cols();
  log.log(1, "Quadratic objective: %zu squares of integer columns (%zu binary), solved as a MILP with tangent rows\n",
          sq.num_terms(), sq.num_terms() - sq.num_general());
  SolveResult r;
  const char* center_env = std::getenv("MILP_SQ_CENTER");  // 0: no relaxation, centers at 0 (in the bounds)
  if (sq.num_general() > 0 && !(center_env && *center_env == '0')) {
    // Centers: the continuous relaxation's point.
    Model rel = model;
    for (VarType& t : rel.col_type) {
      if (t == VarType::Integer || t == VarType::Binary) t = VarType::Continuous;
    }
    const Logger quiet(log.verbosity() >= 2 ? 1 : 0);
    const SolveResult rr = solve_quadratic(rel, params, {}, timer, quiet, nullptr);
    if (rr.has_solution) sq.center_at(rr.x);
    log.log(1, "  relaxation %s %.10g (%.2fs); %zu tangent rows\n", to_string(rr.status), rr.objective,
            timer.elapsed(), sq.num_tangents());
  }
  const char* rounds_env = std::getenv("MILP_SQ_ROUNDS");
  const int max_rounds = rounds_env ? std::max(1, std::atoi(rounds_env)) : 10;
  std::vector<double> best;
  double best_obj = kInf;      // minimization sense
  double best_bound = -kInf;   // minimization sense, valid for the model
  bool converged = false, infeasible = false;
  Status last = Status::Unknown;
  long long nodes = 0, iters = 0;
  for (int round = 0; round < max_rounds; ++round) {
    const Model milp = sq.build();
    std::vector<double> st;
    if (!best.empty()) {
      st = sq.extend(best);
    } else if (!start.empty() && static_cast<Index>(start.size()) == n) {
      bool complete = true;
      for (double v : start) complete = complete && std::isfinite(v);
      st = complete ? sq.extend(start) : start;
      st.resize(static_cast<size_t>(milp.num_cols()), std::numeric_limits<double>::quiet_NaN());
    }
    const SolveResult s = solve_plain(milp, {}, params, st, timer, log, nullptr, tracer, false);
    nodes += s.stats.nodes;
    iters += s.stats.lp_iterations;
    last = s.status;
    if (s.status == Status::Optimal || s.status == Status::TimeLimit || s.status == Status::NodeLimit ||
        s.status == Status::Interrupted || s.status == Status::IterationLimit || s.status == Status::MemoryLimit) {
      best_bound = std::max(best_bound, sense * s.bound);
    }
    if (s.status == Status::Infeasible) {
      infeasible = best.empty();
      converged = !best.empty();
      break;
    }
    if (!s.has_solution) break;
    std::vector<double> x(s.x.begin(), s.x.begin() + n);
    const double obj = sense * model.objective_value(x);
    if (obj < best_obj) best_obj = obj, best = std::move(x);
    if (s.status != Status::Optimal) break;
    const int missed = sq.refine(std::vector<double>(s.x.begin(), s.x.begin() + n));
    if (missed == 0) {
      converged = true;
      break;
    }
    log.log(1, "Tangent round %d: %d squares outside their tangents (MILP %.10g, exact %.10g); %zu tangent rows\n",
            round + 1, missed, s.objective, sense * obj, sq.num_tangents());
  }
  r.stats.nodes = nodes;
  r.stats.lp_iterations = iters;
  if (infeasible) {
    r.status = Status::Infeasible;
    r.bound = sense * kInf;
  } else {
    if (!best.empty()) {
      r.has_solution = true;
      r.x = std::move(best);
      r.objective = model.objective_value(r.x);
      best_bound = std::min(best_bound, best_obj);
    }
    const double tol = std::max(params.mip_abs_gap, params.mip_rel_gap * std::abs(best_obj));
    if (r.has_solution && (converged || best_obj - best_bound <= tol)) {
      r.status = Status::Optimal;
    } else {
      r.status = last == Status::Optimal ? Status::Unknown : last;
    }
    r.bound = sense * best_bound;
  }
  r.stats.root_bound = r.bound;
  r.stats.time = timer.elapsed();
  return r;
}

// A model with a quadratic objective: the QP presolve (qp/qp_presolve.hpp), the convexity test,
// then the barrier (continuous QP) or the MIQP branch-and-bound. A separable convex objective on integer
// columns goes to solve_square_milp (MILP_SQ_MILP=0: the MIQP search).
SolveResult solve_quadratic(const Model& model, const SolverParams& params, const std::vector<double>& start,
                            const Timer& timer, const Logger& log, IncumbentTracer* tracer) {
  {
    const char* on = std::getenv("MILP_SQ_MILP");
    const char* w = std::getenv("MILP_SQ_WINDOW");
    SquareTerms sq;
    if (!(on && *on == '0') && sq.init(model, w ? std::atoi(w) : 25)) {
      return solve_square_milp(model, sq, params, start, timer, log, tracer);
    }
  }
  SolveResult r;
  MipProblem p;
  if (!MipProblem::from_model(model, &p)) {
    r.status = Status::Infeasible;
    r.stats.time = timer.elapsed();
    return r;
  }
  const double sense = model.sense == ObjSense::Maximize ? -1.0 : 1.0;
  log.log(1, "Quadratic objective: %d nonzeros in the lower triangle of Q\n", p.Q.nnz());
  QpPresolver pre;
  const MipProblem* sp = &p;  // the problem solved: p or the presolved one
  if (params.presolve && p.implications.empty()) {
    const double t0 = timer.elapsed();
    const QpPresolveStatus qps = pre.run(p, log);
    PresolveSummary& ps = r.stats.presolve;
    ps.ran = true;
    ps.rows_before = p.num_rows(), ps.cols_before = p.num_cols(), ps.nnz_before = p.A.nnz();
    if (qps != QpPresolveStatus::Infeasible) {
      ps.rows = pre.reduced().num_rows(), ps.cols = pre.reduced().num_cols(), ps.nnz = pre.reduced().A.nnz();
    }
    ps.time = timer.elapsed() - t0;
    if (qps == QpPresolveStatus::Infeasible) {
      r.status = Status::Infeasible;
      r.bound = sense * kInf;
      r.stats.root_bound = r.bound;
      r.stats.time = timer.elapsed();
      return r;
    }
    sp = &pre.reduced();
  }
  // Convexity of the problem solved (cross terms with fixed columns are linear terms).
  const ConvexityResult cv = check_convexity(sp->Q);
  if (!cv.convex) {
    std::string msg = cv.message;
    if (model.sense == ObjSense::Maximize) msg += " (maximization needs a concave objective)";
    log.log(1, "ERROR: %s\n", msg.c_str());
    r.status = Status::NonConvex;
    r.stats.time = timer.elapsed();
    return r;
  }
  if (cv.perturbed) log.log(1, "WARNING: %s\n", cv.message.c_str());
  const std::vector<double>& shift = cv.diag_shift;
  auto postsolve = [&](const std::vector<double>& x) { return sp == &p ? x : pre.postsolve(x); };
  const int threads = params.threads > 0 ? params.threads : static_cast<int>(std::thread::hardware_concurrency());
  if (model.num_integers() == 0) {
    if (!start.empty()) log.log(1, "QP: the start solution is not used (the barrier has no warm start)\n");
    QpSolution qs;
    if (sp->num_cols() == 0) {
      // Presolve fixed every column (and removed every row).
      qs.status = QpStatus::Optimal;
    } else {
      QpOptions qo;
      qo.time_limit = params.time_limit - timer.elapsed();
      qo.threads = std::max(1, threads);
      qo.log = &log;
      qs = solve_qp(*sp, sp->col_lb, sp->col_ub, qo, shift);
    }
    r.stats.barrier_iterations = qs.iterations;
    r.stats.lp_iterations = qs.iterations;
    switch (qs.status) {
      case QpStatus::Optimal: r.status = Status::Optimal; break;
      case QpStatus::Infeasible: r.status = Status::Infeasible; break;
      case QpStatus::Unbounded: r.status = Status::Unbounded; break;
      case QpStatus::InfeasibleOrUnbounded: r.status = Status::InfeasibleOrUnbounded; break;
      case QpStatus::TimeLimit: r.status = Status::TimeLimit; break;
      case QpStatus::NumericalError: r.status = Status::NumericalError; break;
    }
    if (qs.status == QpStatus::Optimal) {
      r.has_solution = true;
      r.x = postsolve(qs.x);
      r.objective = model.objective_value(r.x);
      r.bound = r.objective;
      if (qs.approximate) log.log(1, "QP: solved only to about 1e-6 (numerical difficulties)\n");
      // Duals in the model's sense (p minimizes sense * f).
      if (sp == &p) {
        r.row_dual = qs.y;
        r.reduced_cost = qp_reduced_costs(p, r.x, r.row_dual);
      } else {
        pre.postsolve_duals(r.x, qs.y, r.row_dual, r.reduced_cost);
      }
      r.has_duals = true;
      for (double& v : r.row_dual) v *= sense;
      for (double& v : r.reduced_cost) v *= sense;
    } else if (qs.status == QpStatus::Infeasible) {
      r.bound = sense * kInf;
    }
    r.stats.root_bound = r.bound;
    r.stats.time = timer.elapsed();
    return r;
  }
  MiqpOptions mo;
  mo.threads = std::max(1, threads);
  mo.diag_shift = shift;
  if (tracer) mo.on_trace = tracer->hook();
  if (!start.empty()) mo.start = sp == &p ? start : pre.crush(start);
  MiqpResult mr;
  if (sp->num_cols() == 0) {
    mr.status = Status::Optimal;
    mr.has_solution = true;
    mr.bound = mr.root_bound = sp->offset;
  } else {
    mr = miqp_solve(*sp, params, log, timer, mo);
  }
  r.status = mr.status;
  r.has_solution = mr.has_solution;
  if (mr.has_solution) {
    r.x = postsolve(mr.x);
    r.objective = model.objective_value(r.x);
  }
  r.bound = sense * mr.bound;
  r.stats.nodes = mr.nodes;
  r.stats.barrier_iterations = mr.barrier_iterations;
  r.stats.lp_iterations = mr.barrier_iterations + mr.lp_iterations;
  r.stats.root_bound = sense * mr.root_bound;
  r.stats.time = timer.elapsed();
  return r;
}

}  // namespace

SolveResult solve(const Model& model, const SolverParams& params) { return solve(model, params, {}); }

namespace {

SolveResult solve_model(const Model& model, const SolverParams& params, const std::vector<double>& start,
                        IncumbentReporter* reporter, IncumbentTracer* tracer) {
  Timer timer;
  Logger log(params.verbosity);
  if (model.has_quadratic()) {
    if (!model.has_logical()) return solve_quadratic(model, params, start, timer, log, tracer);
    const LogicalReformulation R = reformulate_logical(model, params.logic_big_m);
    if (!R.implications.empty()) {
      throw std::runtime_error("logical constraints without finite big-M bounds are not supported yet with a "
                               "quadratic objective");
    }
    log.log(1, "Logical constraints: %d semi-continuous, %zu SOS, %zu indicator -> %d auxiliary columns, %d rows\n",
            model.num_semi(), model.sos.size(), model.indicators.size(), R.aux_cols, R.aux_rows);
    SolveResult r = solve_quadratic(R.model, params, start.empty() ? start : R.extend(start), timer, log, tracer);
    if (r.has_solution) {
      r.x.resize(static_cast<size_t>(model.num_cols()));
      snap_logical(model, r.x);
      r.objective = model.objective_value(r.x);
    }
    return r;
  }
  if (!model.has_logical()) return solve_plain(model, {}, params, start, timer, log, reporter, tracer, true);
  // Semi-continuous columns, SOS sets, indicator constraints: solved as a plain MILP with auxiliary
  // binaries (reformulate.hpp). MILP_LOGIC_BIGM overrides params.logic_big_m (diagnostics).
  const char* env_big_m = std::getenv("MILP_LOGIC_BIGM");
  const LogicalReformulation R = reformulate_logical(model, env_big_m ? std::atof(env_big_m) : params.logic_big_m);
  log.log(1, "Logical constraints: %d semi-continuous, %zu SOS, %zu indicator -> %d auxiliary columns, %d rows (%d with big-M), "
          "%zu implications\n", model.num_semi(), model.sos.size(), model.indicators.size(), R.aux_cols, R.aux_rows,
          R.big_m, R.implications.size());
  const std::vector<double> st = start.empty() ? start : R.extend(start);
  SolveResult r = solve_plain(R.model, R.implications, params, st, timer, log, reporter, tracer, false);
  if (r.has_solution) {
    r.x.resize(static_cast<size_t>(model.num_cols()));
    snap_logical(model, r.x);
    r.objective = model.objective_value(r.x);
    const double viol = model.violation(r.x).max();
    if (!(viol <= 1e-5)) log.log(1, "WARNING: solution violates the logical constraints by %.3g\n", viol);
  }
  r.stats.time = timer.elapsed();
  return r;
}

}  // namespace

SolveResult solve(const Model& model, const SolverParams& params_in, const std::vector<double>& start) {
  // The interrupt flags of this solve: the caller's (never written here) and a local one, set when
  // the incumbent callback asks to stop (so the caller's flag stays clear for its next solve).
  std::atomic<bool> local_flag{false};
  SolverParams params = params_in;
  params.interrupt = &local_flag;
  params.on_incumbent = nullptr;  // reported through the reporter only (not copied into sub-solves)
  InterruptScope scope(&local_flag, params_in.interrupt);
  SolverParams cb_params = params;
  cb_params.on_incumbent = params_in.on_incumbent;
  std::optional<IncumbentReporter> reporter;
  if (params_in.on_incumbent && !model.has_quadratic()) reporter.emplace(model, cb_params);
  std::optional<IncumbentTracer> tracer;
  if (params_in.trace_incumbents) tracer.emplace(model.sense);
  SolveResult r = solve_model(model, params, start, reporter ? &*reporter : nullptr, tracer ? &*tracer : nullptr);
  if (tracer) {
    if (r.has_solution) tracer->record(r.objective, "final", r.stats.time);
    r.incumbent_trace = tracer->take();
  }
  // The final solution (postsolved and polished) if it improves on the last one reported.
  if (reporter && r.has_solution) reporter->report(r.x, model.sense == ObjSense::Maximize ? -r.bound : r.bound);
  if (local_flag.load() || (params_in.interrupt && params_in.interrupt->load())) {
    switch (r.status) {
      case Status::TimeLimit:
      case Status::IterationLimit:
      case Status::Unknown:
      case Status::NumericalError:
        r.status = Status::Interrupted;
        break;
      default:
        break;
    }
  }
  return r;
}

}  // namespace milp
