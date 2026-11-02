#pragma once

#include <chrono>

namespace milp {

class Timer {
 public:
  Timer() : start_(Clock::now()) {}
  void reset() { start_ = Clock::now(); }
  double elapsed() const { return std::chrono::duration<double>(Clock::now() - start_).count(); }

 private:
  using Clock = std::chrono::steady_clock;
  Clock::time_point start_;
};

}  // namespace milp
