// Test runner: `popso_tests [substring]` runs all tests whose name contains
// the substring. Returns nonzero if any test fails.

#include <cstring>
#include <exception>

#include "test_framework.hpp"

int main(int argc, char** argv) {
  const char* filter = argc > 1 ? argv[1] : "";
  int run = 0, failed = 0;
  for (const auto& t : milp_test::registry()) {
    if (*filter && !std::strstr(t.name, filter)) continue;
    ++run;
    const int before = milp_test::context().failures;
    try {
      t.fn();
    } catch (const milp_test::RequireFailed&) {
      // already reported
    } catch (const std::exception& e) {
      milp_test::report_failure(t.file, 0, std::string("unexpected exception: ") + e.what());
    } catch (...) {
      milp_test::report_failure(t.file, 0, "unexpected unknown exception");
    }
    const bool ok = milp_test::context().failures == before;
    if (!ok) ++failed;
    std::printf("[%s] %s\n", ok ? " OK " : "FAIL", t.name);
  }
  std::printf("\n%d tests, %d passed, %d failed (%d checks, %d failed)\n", run, run - failed, failed,
              milp_test::context().checks, milp_test::context().failures);
  if (run == 0) std::printf("no tests matched filter '%s'\n", filter);
  return failed == 0 && run > 0 ? 0 : 1;
}
