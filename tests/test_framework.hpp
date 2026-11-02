#pragma once
// Tiny std-only unit test harness. Tests register themselves at static
// initialization time:
//
//   TEST(sparse_transpose) { CHECK(x == 1); CHECK_NEAR(a, b, 1e-9); REQUIRE(ok); }
//
// CHECK* record a failure and continue; REQUIRE aborts the current test.

#include <cmath>
#include <cstdio>
#include <exception>
#include <functional>
#include <string>
#include <vector>

namespace milp_test {

struct TestCase {
  const char* name;
  const char* file;
  void (*fn)();
};

inline std::vector<TestCase>& registry() {
  static std::vector<TestCase> tests;
  return tests;
}

struct Registrar {
  Registrar(const char* name, const char* file, void (*fn)()) { registry().push_back({name, file, fn}); }
};

struct Context {
  int checks = 0;
  int failures = 0;
};

inline Context& context() {
  static Context ctx;
  return ctx;
}

struct RequireFailed {};

inline void report_failure(const char* file, int line, const std::string& what) {
  ++context().failures;
  std::fprintf(stderr, "  %s:%d: FAILED: %s\n", file, line, what.c_str());
}

inline bool check(bool ok, const char* file, int line, const char* expr) {
  ++context().checks;
  if (!ok) report_failure(file, line, expr);
  return ok;
}

inline bool check_near(double a, double b, double tol, const char* file, int line, const char* ea, const char* eb) {
  ++context().checks;
  const bool ok = (a == b) || std::abs(a - b) <= tol;
  if (!ok) {
    char buf[256];
    std::snprintf(buf, sizeof(buf), "%s ~= %s (%.17g vs %.17g, tol %g)", ea, eb, a, b, tol);
    report_failure(file, line, buf);
  }
  return ok;
}

}  // namespace milp_test

#define MILP_TEST_CAT2(a, b) a##b
#define MILP_TEST_CAT(a, b) MILP_TEST_CAT2(a, b)

#define TEST(name)                                                                            \
  static void MILP_TEST_CAT(milp_test_fn_, name)();                                           \
  static const ::milp_test::Registrar MILP_TEST_CAT(milp_test_reg_, name)(#name, __FILE__,    \
                                                                          &MILP_TEST_CAT(milp_test_fn_, name)); \
  static void MILP_TEST_CAT(milp_test_fn_, name)()

#define CHECK(cond) ::milp_test::check(static_cast<bool>(cond), __FILE__, __LINE__, #cond)

#define CHECK_EQ(a, b) ::milp_test::check((a) == (b), __FILE__, __LINE__, #a " == " #b)

#define CHECK_NEAR(a, b, tol) \
  ::milp_test::check_near(static_cast<double>(a), static_cast<double>(b), (tol), __FILE__, __LINE__, #a, #b)

#define REQUIRE(cond)                                                                    \
  do {                                                                                   \
    if (!::milp_test::check(static_cast<bool>(cond), __FILE__, __LINE__, #cond)) throw ::milp_test::RequireFailed{}; \
  } while (0)

// Passes if `expr` throws an exception of type `ex` (or derived).
#define CHECK_THROWS(expr, ex)                                                            \
  do {                                                                                    \
    bool milp_test_threw = false;                                                         \
    try {                                                                                 \
      (void)(expr);                                                                       \
    } catch (const ex&) {                                                                 \
      milp_test_threw = true;                                                             \
    } catch (...) {                                                                       \
    }                                                                                     \
    ::milp_test::check(milp_test_threw, __FILE__, __LINE__, "throws " #ex ": " #expr);   \
  } while (0)
