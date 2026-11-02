#pragma once
// Minimal thread-safe logger writing to stdout.

#include <cstdarg>
#include <cstdio>
#include <mutex>

namespace milp {

class Logger {
 public:
  explicit Logger(int verbosity = 1) : verbosity_(verbosity) {}
  int verbosity() const { return verbosity_; }
  void set_verbosity(int v) { verbosity_ = v; }

#if defined(__GNUC__) || defined(__clang__)
  __attribute__((format(printf, 3, 4)))
#endif
  void log(int level, const char* fmt, ...) const {
    if (level > verbosity_) return;
    std::lock_guard<std::mutex> lock(mutex_);
    va_list args;
    va_start(args, fmt);
    std::vfprintf(stdout, fmt, args);
    va_end(args);
    std::fflush(stdout);
  }

 private:
  int verbosity_;
  mutable std::mutex mutex_;
};

}  // namespace milp
