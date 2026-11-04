#pragma once
// gzip decompression (RFC 1952 container, RFC 1951 deflate), standard library only.
//
// The decoder pulls compressed input from a source in chunks and pushes decompressed data to a
// sink in chunks, keeping only the 32 KB deflate window plus one output chunk in memory. Every
// member's CRC-32 and length are checked; concatenated members are decoded one after another.
// Corrupt or truncated data throws GzipError.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>

namespace milp {

class GzipError : public std::runtime_error {
 public:
  using std::runtime_error::runtime_error;
};

// CRC-32 (IEEE 802.3, as in gzip) of data, continuing from crc (0 to start).
uint32_t crc32_update(uint32_t crc, const unsigned char* data, size_t n);

// True if the data starts with the gzip magic bytes 1f 8b.
bool is_gzip(const void* data, size_t n);

// Source: fills buf with up to cap bytes, returns how many (0 at the end of the input).
using GzipSource = std::function<size_t(unsigned char* buf, size_t cap)>;
// Sink: receives the decompressed data in order.
using GzipSink = std::function<void(const char* data, size_t n)>;

// Decompresses a whole gzip stream (one or more members). Returns the number of bytes produced.
uint64_t gunzip(const GzipSource& source, const GzipSink& sink);

// Cap on the decompressed size kept in memory by the functions below (deflate reaches ratios
// of about 1000:1, so a small corrupt or malicious file could otherwise exhaust memory). Process-
// wide; 0 restores the default, 16 GiB.
constexpr uint64_t kDefaultGzipMaxOutput = uint64_t{16} << 30;
void set_gzip_max_output(uint64_t bytes);
uint64_t gzip_max_output();

// Convenience: decompresses a gzip stream held in memory. Throws GzipError once the output would
// exceed max_output bytes.
std::string gunzip(const std::string& compressed, uint64_t max_output = gzip_max_output());

// Reads a file into memory, decompressing it if it is gzip data (detected by the magic bytes,
// whatever the file name). Throws std::runtime_error (GzipError for corrupt data, or once the
// decompressed data would exceed max_output bytes) with the path.
std::string read_maybe_gzipped_file(const std::string& path, uint64_t max_output = gzip_max_output());

}  // namespace milp
