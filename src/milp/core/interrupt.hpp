#pragma once
// User interrupts of a running solve.
//
// solve() (solver/solver.hpp) installs the solve's interrupt flags (SolverParams::interrupt and its
// own, set by the incumbent callback) as the calling thread's current flags; every thread the
// solver starts inherits them (wrap the thread body in inherit_interrupt). Long-running
// loops (simplex and barrier iterations, the branch-and-bound loop, presolve rounds, heuristics)
// poll interrupted() and stop as if a time limit had been reached; solve() then reports
// Status::Interrupted. Setting the flag is a single atomic store, so it may be done from any
// thread or a signal handler (std::atomic<bool> is lock-free on every supported platform).
//
// The flag is per solve (thread-local pointer), so separate models solved concurrently in
// different threads do not interrupt each other.
//
// A third flag belongs to the running branch-and-bound search (SearchStopScope): the search sets it
// once it has solved the problem (an incumbent meets the trivial objective bound), which stops the
// LPs, cut loop and heuristics still running for it the same way.

#include <atomic>
#include <chrono>

namespace milp {

namespace detail {
// The current thread's stop state: up to two interrupt flags (solve() combines the caller's flag
// with its own, which the incumbent callback sets), the search's flag (see SearchStopScope) and a
// deadline (see DeadlineScope).
struct InterruptState {
  const std::atomic<bool>* flag = nullptr;
  const std::atomic<bool>* flag2 = nullptr;
  const std::atomic<bool>* search = nullptr;
  const std::atomic<bool>* local = nullptr;  // see LocalStopScope
  std::chrono::steady_clock::time_point deadline = std::chrono::steady_clock::time_point::max();
};
inline thread_local InterruptState t_interrupt;
}  // namespace detail

// True once the current solve has been interrupted (cheap: thread-local loads and relaxed atomic
// loads). Ignores the deadline.
inline bool interrupted() {
  const detail::InterruptState& s = detail::t_interrupt;
  return (s.flag != nullptr && s.flag->load(std::memory_order_relaxed)) ||
         (s.flag2 != nullptr && s.flag2->load(std::memory_order_relaxed)) ||
         (s.search != nullptr && s.search->load(std::memory_order_relaxed)) ||
         (s.local != nullptr && s.local->load(std::memory_order_relaxed));
}

// True once the current deadline (DeadlineScope) has passed; reads the clock only when a deadline
// is set.
inline bool deadline_passed() {
  const auto d = detail::t_interrupt.deadline;
  return d != std::chrono::steady_clock::time_point::max() && std::chrono::steady_clock::now() >= d;
}

// interrupted() or deadline_passed(): for loops that poll interrupts at a fine grain (inside a
// factorization) and should honour a time limit just as fast.
inline bool interrupted_or_expired() { return interrupted() || deadline_passed(); }

// Makes `flag` (and `flag2`, if any) the current thread's interrupt flags for the scope's
// lifetime. The search flag and the deadline are kept.
class InterruptScope {
 public:
  explicit InterruptScope(const std::atomic<bool>* flag, const std::atomic<bool>* flag2 = nullptr)
      : saved_(detail::t_interrupt) {
    detail::t_interrupt.flag = flag;
    detail::t_interrupt.flag2 = flag2;
  }
  ~InterruptScope() { detail::t_interrupt = saved_; }
  InterruptScope(const InterruptScope&) = delete;
  InterruptScope& operator=(const InterruptScope&) = delete;

 private:
  detail::InterruptState saved_;
};

// Makes `flag` the current thread's search flag for the scope's lifetime, unless one is set already
// (a search inside another one, e.g. a sub-MIP, keeps the outer search's flag, which stops it too).
// The other flags and the deadline are kept.
class SearchStopScope {
 public:
  explicit SearchStopScope(const std::atomic<bool>* flag) : saved_(detail::t_interrupt.search) {
    if (saved_ == nullptr) detail::t_interrupt.search = flag;
  }
  ~SearchStopScope() { detail::t_interrupt.search = saved_; }
  SearchStopScope(const SearchStopScope&) = delete;
  SearchStopScope& operator=(const SearchStopScope&) = delete;

 private:
  const std::atomic<bool>* saved_;
};

// Makes `flag` the current thread's local stop flag for the scope's lifetime (also inside a search
// that has a search flag already): stops work a search runs beside its tree (the root heuristics in
// the background) once the tree has ended. The other flags and the deadline are kept.
class LocalStopScope {
 public:
  explicit LocalStopScope(const std::atomic<bool>* flag) : saved_(detail::t_interrupt.local) {
    detail::t_interrupt.local = flag;
  }
  ~LocalStopScope() { detail::t_interrupt.local = saved_; }
  LocalStopScope(const LocalStopScope&) = delete;
  LocalStopScope& operator=(const LocalStopScope&) = delete;

 private:
  const std::atomic<bool>* saved_;
};

// Sets the current thread's deadline to `seconds` from now (no later than an enclosing one; an
// infinite or NaN value keeps the enclosing deadline) for the scope's lifetime. Threads started
// with inherit_interrupt inside the scope inherit it.
class DeadlineScope {
 public:
  explicit DeadlineScope(double seconds) : saved_(detail::t_interrupt.deadline) {
    using clock = std::chrono::steady_clock;
    if (!(seconds < 1e9)) return;  // ~30 years: no deadline (also NaN)
    const auto d = clock::now() + std::chrono::duration_cast<clock::duration>(
                                      std::chrono::duration<double>(seconds > 0.0 ? seconds : 0.0));
    if (d < saved_) detail::t_interrupt.deadline = d;
  }
  ~DeadlineScope() { detail::t_interrupt.deadline = saved_; }
  DeadlineScope(const DeadlineScope&) = delete;
  DeadlineScope& operator=(const DeadlineScope&) = delete;

 private:
  std::chrono::steady_clock::time_point saved_;
};

// Wraps a thread body so that it runs with the creating thread's interrupt flags and deadline:
// std::thread t(inherit_interrupt([&] { ... }));
template <class F>
auto inherit_interrupt(F f) {
  return [state = detail::t_interrupt, f = std::move(f)]() mutable {
    const detail::InterruptState saved = detail::t_interrupt;
    detail::t_interrupt = state;
    f();
    detail::t_interrupt = saved;
  };
}

}  // namespace milp
