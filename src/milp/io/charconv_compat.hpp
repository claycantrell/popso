#pragma once
// Fast, locale-independent number <-> text conversion with a fallback for
// standard libraries that lack floating-point <charconv> (e.g. libc++ < 20
// has no floating-point from_chars, libc++ < 14 no floating-point to_chars).

#include <charconv>
#include <cstdio>
#include <cstdlib>
#include <string_view>
#include <system_error>
#include <version>

#if defined(MILP_NO_FP_CHARCONV)  // force the fallback (for testing)
#define MILP_HAS_FP_FROM_CHARS 0
#define MILP_HAS_FP_TO_CHARS 0
#elif defined(__cpp_lib_to_chars) && __cpp_lib_to_chars >= 201611L
#define MILP_HAS_FP_FROM_CHARS 1
#define MILP_HAS_FP_TO_CHARS 1
#elif defined(_LIBCPP_VERSION)
#define MILP_HAS_FP_FROM_CHARS (_LIBCPP_VERSION >= 200000)
#define MILP_HAS_FP_TO_CHARS (_LIBCPP_VERSION >= 14000)
#elif defined(__GLIBCXX__) && defined(_GLIBCXX_RELEASE) && _GLIBCXX_RELEASE >= 11
#define MILP_HAS_FP_FROM_CHARS 1
#define MILP_HAS_FP_TO_CHARS 1
#elif defined(_MSC_VER)
#define MILP_HAS_FP_FROM_CHARS 1
#define MILP_HAS_FP_TO_CHARS 1
#else
#define MILP_HAS_FP_FROM_CHARS 0
#define MILP_HAS_FP_TO_CHARS 0
#endif

namespace milp::io_detail {

// Parses the whole token as a double. Accepts a leading '+', decimal and
// exponent notation and "inf"/"infinity"/"nan" (any case). Returns false if
// the token is empty or not entirely consumed.
//
// The fallback path uses strtod, which requires that the character after the
// token is not part of a number: callers pass tokens from a buffer where every
// token is followed by whitespace or a terminating NUL.
inline bool parse_double(std::string_view tok, double& out) {
  if (!tok.empty() && tok.front() == '+') tok.remove_prefix(1);
  if (tok.empty()) return false;
#if MILP_HAS_FP_FROM_CHARS
  auto res = std::from_chars(tok.data(), tok.data() + tok.size(), out);
  return res.ec == std::errc() && res.ptr == tok.data() + tok.size();
#else
  if (tok.front() == ' ' || tok.front() == '\t') return false;
  char* end = nullptr;
  out = std::strtod(tok.data(), &end);
  return end == tok.data() + tok.size();
#endif
}

// Writes the shortest representation that round-trips exactly. `buf` must
// hold at least 32 chars. Returns the number of chars written.
inline int format_double(char* buf, double v) {
#if MILP_HAS_FP_TO_CHARS
  auto res = std::to_chars(buf, buf + 32, v);
  return static_cast<int>(res.ptr - buf);
#else
  return std::snprintf(buf, 32, "%.17g", v);
#endif
}

}  // namespace milp::io_detail
