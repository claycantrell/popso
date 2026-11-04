// Interrupting a solve (SolverParams::interrupt) and the incumbent callback (SolverParams::on_incumbent).

#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <thread>
#include <vector>

#include "milp/core/interrupt.hpp"
#include "milp/core/model.hpp"
#include "milp/io/mps_reader.hpp"
#include "milp/solver/solver.hpp"
#include "test_framework.hpp"

using namespace milp;

namespace {

bool load(const std::string& name, Model* m) {
  try {
    *m = read_mps(std::string(MILP_SOURCE_DIR) + "/benchmarks/instances/" + name + ".mps");
    return true;
  } catch (const std::exception&) {
    return false;  // generated instances are optional
  }
}

SolverParams params2(bool deterministic) {
  SolverParams sp;
  sp.verbosity = 0;
  sp.threads = 2;
  sp.deterministic = deterministic;
  return sp;
}

double seconds_since(std::chrono::steady_clock::time_point t0) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// Latency limits below are for optimized builds; address / undefined-behavior sanitizer builds run
// several times slower (the linux-sanitizers CI job measured 2.2 s for a 0.2 s stop).
#if defined(__SANITIZE_ADDRESS__)
constexpr double kSlow = 4.0;
#elif defined(__has_feature)
#if __has_feature(address_sanitizer)
constexpr double kSlow = 4.0;
#else
constexpr double kSlow = 1.0;
#endif
#else
constexpr double kSlow = 1.0;
#endif

}  // namespace

TEST(interrupt_scope_inheritance) {
  std::atomic<bool> flag{false};
  CHECK(!interrupted());
  {
    InterruptScope scope(&flag);
    CHECK(!interrupted());
    flag = true;
    CHECK(interrupted());
    bool in_thread = false;
    std::thread t(inherit_interrupt([&] { in_thread = interrupted(); }));
    t.join();
    CHECK(in_thread);
    bool plain_thread = true;
    std::thread u([&] { plain_thread = interrupted(); });
    u.join();
    CHECK(!plain_thread);
  }
  CHECK(!interrupted());
}

TEST(interrupt_stops_the_search) {
  Model m;
  if (!load("cflp_hard_1", &m)) return;
  for (int det = 0; det < 2; ++det) {
    // Set before the solve: stops at once.
    std::atomic<bool> flag{true};
    SolverParams sp = params2(det != 0);
    sp.interrupt = &flag;
    auto t0 = std::chrono::steady_clock::now();
    SolveResult r = solve(m, sp);
    CHECK(r.status == Status::Interrupted);
    CHECK(seconds_since(t0) < 1.0 * kSlow);
    // Set from another thread during the search: stops within a fraction of a second, with the
    // incumbent found so far and a valid bound.
    flag = false;
    t0 = std::chrono::steady_clock::now();
    std::thread setter([&] {
      std::this_thread::sleep_for(std::chrono::milliseconds(700));
      flag = true;
    });
    r = solve(m, sp);
    const double t = seconds_since(t0);
    setter.join();
    CHECK(r.status == Status::Interrupted);
    // Typically 3-40 ms after the flag; slow shared CI machines need more room.
    if (!(t < 0.7 + 2.0 * kSlow)) std::fprintf(stderr, "interrupt latency: %.3f s after the flag\n", t - 0.7);
    CHECK(t < 0.7 + 2.0 * kSlow);
    if (r.has_solution) {
      CHECK(m.violation(r.x).max() <= 1e-6);
      CHECK(r.bound <= r.objective + 1e-6);
    }
  }
}

TEST(incumbent_callback_order_and_stop) {
  for (const char* name : {"lotsize_small_1", "gap_small_1", "cflp_small_1"}) {
    Model m;
    if (!load(name, &m)) continue;
    const double sense = m.sense == ObjSense::Maximize ? -1.0 : 1.0;
    std::vector<std::vector<double>> runs;
    for (int rep = 0; rep < 2; ++rep) {
      SolverParams sp = params2(true);
      std::vector<double> objs;
      bool ok = true;
      sp.on_incumbent = [&](const std::vector<double>& x, double obj, double bound) {
        ok = ok && x.size() == static_cast<size_t>(m.num_cols());
        ok = ok && std::abs(m.objective_value(x) - obj) <= 1e-9 * (1.0 + std::abs(obj));
        ok = ok && sense * bound <= sense * obj + 1e-9;
        ok = ok && m.violation(x).max() <= 1e-5;
        if (!objs.empty()) ok = ok && sense * obj < sense * objs.back();  // strictly improving
        objs.push_back(obj);
        return false;
      };
      const SolveResult r = solve(m, sp);
      CHECK(r.status == Status::Optimal);
      CHECK(ok);
      REQUIRE(!objs.empty());
      CHECK_NEAR(objs.back(), r.objective, 1e-6 * (1.0 + std::abs(r.objective)));  // the final solution is reported
      runs.push_back(objs);
    }
    CHECK(runs[0] == runs[1]);  // deterministic mode: the same sequence
    // Asking to stop at the first incumbent: no further calls, status Interrupted (unless the
    // search had finished anyway).
    SolverParams sp = params2(false);
    int calls = 0;
    sp.on_incumbent = [&](const std::vector<double>&, double, double) {
      ++calls;
      return true;
    };
    const SolveResult r = solve(m, sp);
    CHECK(calls == 1);
    CHECK(r.status == Status::Interrupted || r.status == Status::Optimal);
    CHECK(r.has_solution);
  }
}

TEST(interrupt_lp_and_qp) {
  // Pure LP (simplex and barrier) and MIQP paths stop too.
  std::atomic<bool> flag{true};
  for (const char* name : {"randlp_medium_1", "miqprand_medium_1"}) {
    Model m;
    if (!load(name, &m)) continue;
    for (int method = 1; method <= 2; ++method) {
      SolverParams sp = params2(false);
      sp.lp_method = static_cast<LpMethod>(method);
      sp.interrupt = &flag;
      const auto t0 = std::chrono::steady_clock::now();
      const SolveResult r = solve(m, sp);
      CHECK(r.status == Status::Interrupted);
      if (r.status != Status::Interrupted) std::printf("  %s method %d: %s\n", name, method, to_string(r.status));
      CHECK(seconds_since(t0) < 1.0 * kSlow);
    }
  }
}

TEST(incumbent_callback_bound_is_valid) {
  // The bound passed with each incumbent must be a valid bound on the optimum (also in the
  // deterministic mode, whose callbacks come at the merge of a round, while the round's children
  // are not in the node pool yet).
  int checked = 0;
  for (const char* name : {"gap_medium_1", "setcover_medium_1", "randmilp_medium_1", "mkp_small_1"}) {
    Model m;
    if (!load(name, &m)) continue;
    const double sense = m.sense == ObjSense::Maximize ? -1.0 : 1.0;
    for (int det = 0; det < 2; ++det) {
      for (unsigned seed = 0; seed < 2; ++seed) {
        SolverParams sp = params2(det != 0);
        sp.seed = seed;
        sp.time_limit = 20;
        std::vector<double> bounds;
        sp.on_incumbent = [&](const std::vector<double>&, double, double bound) {
          bounds.push_back(bound);
          return false;
        };
        const SolveResult r = solve(m, sp);
        if (r.status != Status::Optimal) continue;
        ++checked;
        for (double b : bounds) {
          const bool ok = sense * b <= sense * r.objective + 1e-6 * (1.0 + std::abs(r.objective));
          if (!ok) std::printf("  %s det %d seed %u: bound %.10g beyond the optimum %.10g\n", name, det, seed, b, r.objective);
          CHECK(ok);
        }
      }
    }
  }
  CHECK(checked > 0);
}

TEST(interrupt_at_random_times_gives_valid_results) {
  // Wherever the interrupt lands (presolve, root LP, cuts, heuristics, tree, duals), the result
  // must stay valid: an Optimal status has the true optimum, a solution is feasible, the bound is
  // valid, and an interrupted solve is never reported infeasible or unbounded.
  const int delays_us[] = {0, 100, 300, 1000, 2000, 4000, 8000, 15000, 30000, 60000, 120000};
  int runs = 0;
  for (const char* name : {"lotsize_small_1", "gap_small_1", "cflp_small_1", "randlp_small_1", "netdesign_small_1",
                           "semilot_small_1", "pwlnet_small_1", "tsp_small_1", "portfolio_small_1",
                           "miqprand_small_1", "transport_small_1", "setcover_small_1"}) {
    Model m;
    if (!load(name, &m)) continue;
    const double sense = m.sense == ObjSense::Maximize ? -1.0 : 1.0;
    SolverParams ref_p = params2(false);
    const SolveResult ref = solve(m, ref_p);
    if (ref.status != Status::Optimal) continue;
    const double tol = 1e-6 * (1.0 + std::abs(ref.objective));
    for (int det = 0; det < 2; ++det) {
      for (int d : delays_us) {
        std::atomic<bool> flag{false};
        SolverParams sp = params2(det != 0);
        sp.duals = 2;
        sp.interrupt = &flag;
        std::thread setter([&] {
          std::this_thread::sleep_for(std::chrono::microseconds(d));
          flag = true;
        });
        const SolveResult r = solve(m, sp);
        setter.join();
        ++runs;
        bool ok = r.status == Status::Optimal || r.status == Status::Interrupted;
        if (r.status == Status::Optimal) ok = ok && std::abs(r.objective - ref.objective) <= 1e-5 * (1.0 + std::abs(ref.objective));
        if (r.has_solution) ok = ok && m.violation(r.x).max() <= 1e-5;
        if (r.has_solution) ok = ok && sense * r.objective >= sense * ref.objective - tol;
        ok = ok && sense * r.bound <= sense * ref.objective + tol;
        if (!ok)
          std::printf("  %s det %d delay %d us: %s obj %.10g bound %.10g (optimum %.10g)\n", name, det, d,
                      to_string(r.status), r.objective, r.bound, ref.objective);
        CHECK(ok);
      }
    }
  }
  CHECK(runs > 0);
}

TEST(interrupt_two_flags_and_deadline) {
  std::atomic<bool> a{false}, b{false};
  {
    InterruptScope scope(&a, &b);
    CHECK(!interrupted());
    b = true;
    CHECK(interrupted());
    bool in_thread = false;
    std::thread t(inherit_interrupt([&] { in_thread = interrupted(); }));
    t.join();
    CHECK(in_thread);
    b = false;
    CHECK(!interrupted());
  }
  // Deadlines: never earlier than an enclosing one is kept, inherited by threads, not an
  // interrupt (interrupted() ignores them; interrupted_or_expired() does not).
  CHECK(!deadline_passed());
  {
    DeadlineScope d(0.0);
    CHECK(deadline_passed());
    CHECK(interrupted_or_expired());
    CHECK(!interrupted());
    {
      DeadlineScope later(100.0);  // the enclosing (earlier) deadline stays
      CHECK(deadline_passed());
    }
    bool in_thread = false;
    std::thread t(inherit_interrupt([&] { in_thread = deadline_passed(); }));
    t.join();
    CHECK(in_thread);
  }
  CHECK(!deadline_passed());
  {
    DeadlineScope none(kInf);
    CHECK(!deadline_passed());
    DeadlineScope far(100.0);
    CHECK(!deadline_passed());
  }
}

TEST(incumbent_callback_stop_leaves_caller_flag_clear) {
  // A stop asked by the incumbent callback must not set the caller's interrupt flag: the same
  // params (flag, callback) are reused for a second solve that must run normally.
  Model m;
  if (!load("cflp_hard_1", &m)) return;
  std::atomic<bool> flag{false};
  SolverParams sp = params2(true);
  sp.interrupt = &flag;
  sp.node_limit = 50;
  bool stop = true;
  int calls = 0;
  sp.on_incumbent = [&](const std::vector<double>&, double, double) {
    ++calls;
    return stop;
  };
  SolveResult r = solve(m, sp);
  CHECK(r.status == Status::Interrupted);
  CHECK(calls == 1);
  CHECK(!flag.load());
  stop = false;
  r = solve(m, sp);
  CHECK(r.status == Status::NodeLimit);
  CHECK(r.stats.nodes >= 50);
  // The caller's flag still interrupts.
  flag = true;
  r = solve(m, sp);
  CHECK(r.status == Status::Interrupted);
}

TEST(barrier_time_limit_is_fine_grained) {
  // The time limit is checked inside the barrier's factorizations (as interrupts are): one IPM
  // iteration of randlp_large takes seconds, the solve must stop much sooner.
  Model m;
  if (!load("randlp_large_1", &m)) return;
  for (int method = 0; method <= 2; method += 2) {
    SolverParams sp = params2(false);
    sp.lp_method = static_cast<LpMethod>(method);
    sp.time_limit = 0.0001;
    const auto t0 = std::chrono::steady_clock::now();
    const SolveResult r = solve(m, sp);
    const double t = seconds_since(t0);
    CHECK(r.status == Status::TimeLimit);
    CHECK(t < 2.0 * kSlow);
    if (t >= 2.0 * kSlow) std::printf("  method %d: %.2fs\n", method, t);
  }
}
