#include "milp/io/gzip.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <memory>
#include <system_error>
#include <vector>

namespace milp {

namespace {

const std::array<uint32_t, 256>& crc_table() {
  static const std::array<uint32_t, 256> table = [] {
    std::array<uint32_t, 256> t{};
    for (uint32_t i = 0; i < 256; ++i) {
      uint32_t c = i;
      for (int k = 0; k < 8; ++k) c = (c & 1) ? 0xEDB88320u ^ (c >> 1) : c >> 1;
      t[i] = c;
    }
    return t;
  }();
  return table;
}

// Canonical Huffman code as a single-level lookup table indexed by the next `bits` input bits
// (deflate sends codes most significant bit first, so the table is indexed by the bit-reversed
// code). Entry: symbol << 4 | code length; length 0 marks bit patterns of an incomplete code.
struct Huffman {
  std::vector<uint16_t> table;
  int bits = 0;

  void build(const uint8_t* lengths, int n) {
    int count[16] = {};
    int max_len = 0;
    for (int s = 0; s < n; ++s) {
      if (lengths[s] > 15) throw GzipError("invalid Huffman code length");
      ++count[lengths[s]];
      max_len = std::max(max_len, static_cast<int>(lengths[s]));
    }
    count[0] = 0;
    // Over-subscribed codes are invalid; incomplete ones are accepted (an unused pattern is an
    // error only when it occurs in the data).
    int left = 1;
    for (int len = 1; len <= 15; ++len) {
      left <<= 1;
      left -= count[len];
      if (left < 0) throw GzipError("over-subscribed Huffman code");
    }
    int next[16] = {};
    int code = 0;
    for (int len = 1; len <= 15; ++len) {
      code = (code + count[len - 1]) << 1;
      next[len] = code;
    }
    bits = std::max(max_len, 1);
    table.assign(static_cast<size_t>(1) << bits, 0);
    for (int s = 0; s < n; ++s) {
      const int len = lengths[s];
      if (len == 0) continue;
      const int c = next[len]++;
      int r = 0;
      for (int k = 0; k < len; ++k) r |= ((c >> k) & 1) << (len - 1 - k);
      const uint16_t entry = static_cast<uint16_t>(s << 4 | len);
      for (size_t k = static_cast<size_t>(r); k < table.size(); k += static_cast<size_t>(1) << len) table[k] = entry;
    }
  }
};

constexpr uint16_t kLenBase[29] = {3,  4,  5,  6,  7,  8,  9,  10, 11,  13,  15,  17,  19,  23, 27,
                                   31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258};
constexpr uint8_t kLenExtra[29] = {0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0};
constexpr uint16_t kDistBase[30] = {1,    2,    3,    4,    5,    7,     9,     13,    17,  25,
                                    33,   49,   65,   97,   129,  193,   257,   385,   513, 769,
                                    1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577};
constexpr uint8_t kDistExtra[30] = {0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13};
constexpr uint8_t kClOrder[19] = {16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15};

class Inflater {
 public:
  Inflater(const GzipSource& source, const GzipSink& sink)
      : source_(source), sink_(sink), in_(kInChunk), out_(new unsigned char[kOutCap]) {
    uint8_t len[288];
    for (int s = 0; s < 288; ++s) len[s] = s < 144 ? 8 : s < 256 ? 9 : s < 280 ? 7 : 8;
    fixed_lit_.build(len, 288);
    uint8_t dl[30];
    std::fill(dl, dl + 30, static_cast<uint8_t>(5));
    fixed_dist_.build(dl, 30);
  }

  uint64_t run() {
    refill();
    if (bitcnt_ < 16) throw GzipError("not gzip data (too short)");
    for (;;) {
      member_ = 0;
      crc_ = 0;
      header();
      for (bool last = false; !last;) {
        last = take(1) != 0;
        const uint32_t type = take(2);
        if (type == 0) stored();
        else if (type == 1) codes(fixed_lit_, fixed_dist_);
        else if (type == 2) dynamic();
        else throw GzipError("invalid deflate block type");
      }
      flush(true);
      trailer();
      // Another member, or the end (zero padding after the last member is tolerated).
      refill();
      if (bitcnt_ == 0) break;
      if (bitcnt_ >= 16 && (bitbuf_ & 0xFFFF) == 0x8B1F) continue;
      while (bitcnt_ > 0) {
        if ((bitbuf_ & 0xFF) != 0) throw GzipError("trailing garbage after the gzip data");
        drop(8);
        refill();
      }
      break;
    }
    return total_;
  }

 private:
  static constexpr size_t kInChunk = 1 << 16;
  static constexpr size_t kWindow = 32768;
  static constexpr size_t kOutChunk = 1 << 18;
  static constexpr size_t kOutCap = kWindow + kOutChunk + 512;

  // --- input bits (LSB first)
  void refill() {
    while (bitcnt_ <= 56) {
      if (in_pos_ == in_len_) {
        if (in_eof_) return;
        in_len_ = source_(in_.data(), in_.size());
        in_pos_ = 0;
        if (in_len_ == 0) {
          in_eof_ = true;
          return;
        }
      }
      bitbuf_ |= static_cast<uint64_t>(in_[in_pos_++]) << bitcnt_;
      bitcnt_ += 8;
    }
  }
  uint32_t take(int n) {
    if (bitcnt_ < n) {
      refill();
      if (bitcnt_ < n) throw GzipError("unexpected end of gzip data");
    }
    const uint32_t v = static_cast<uint32_t>(bitbuf_ & ((uint64_t{1} << n) - 1));
    drop(n);
    return v;
  }
  void drop(int n) {
    bitbuf_ >>= n;
    bitcnt_ -= n;
  }
  void align() { drop(bitcnt_ & 7); }
  uint32_t byte() { return take(8); }
  uint32_t le16() {
    const uint32_t lo = byte();
    return lo | byte() << 8;
  }
  uint32_t le32() {
    const uint32_t lo = le16();
    return lo | le16() << 16;
  }

  int decode(const Huffman& h) {
    if (bitcnt_ < h.bits) refill();
    const uint16_t e = h.table[static_cast<size_t>(bitbuf_ & ((uint64_t{1} << h.bits) - 1))];
    const int len = e & 15;
    if (len == 0) throw GzipError("invalid Huffman code in deflate data");
    if (len > bitcnt_) throw GzipError("unexpected end of gzip data");
    drop(len);
    return e >> 4;
  }

  // --- gzip member header and trailer (RFC 1952)
  void header() {
    if (byte() != 0x1F || byte() != 0x8B) throw GzipError("not gzip data (bad magic bytes)");
    if (byte() != 8) throw GzipError("unsupported gzip compression method");
    const uint32_t flg = byte();
    if (flg & 0xE0) throw GzipError("reserved gzip header flags set");
    for (int k = 0; k < 6; ++k) byte();  // MTIME, XFL, OS
    if (flg & 4) {                        // FEXTRA
      const uint32_t xlen = le16();
      for (uint32_t k = 0; k < xlen; ++k) byte();
    }
    if (flg & 8) {  // FNAME
      while (byte() != 0) {
      }
    }
    if (flg & 16) {  // FCOMMENT
      while (byte() != 0) {
      }
    }
    if (flg & 2) le16();  // FHCRC (not checked)
  }
  void trailer() {
    align();
    const uint32_t crc = le32(), isize = le32();
    if (crc != crc_) throw GzipError("gzip CRC-32 mismatch (corrupt data)");
    if (isize != static_cast<uint32_t>(member_)) throw GzipError("gzip length mismatch (corrupt data)");
  }

  // --- output with the sliding window
  void flush(bool all) {
    const size_t n = pos_ - flushed_;
    if (n > 0) {
      crc_ = crc32_update(crc_, out_.get() + flushed_, n);
      sink_(reinterpret_cast<const char*>(out_.get() + flushed_), n);
      flushed_ = pos_;
    }
    if (all || pos_ > kWindow) {
      // Keep the last 32 KB for back references.
      const size_t keep = std::min(pos_, kWindow);
      std::memmove(out_.get(), out_.get() + pos_ - keep, keep);
      pos_ = flushed_ = keep;
    }
    if (all) pos_ = flushed_ = 0;  // a new member starts without history
  }
  void room() {
    if (pos_ + 258 > kOutCap) flush(false);
  }
  void put(unsigned char c) {
    out_[pos_++] = c;
    ++member_;
    ++total_;
  }

  // --- deflate blocks (RFC 1951)
  void stored() {
    align();
    const uint32_t len = le16(), nlen = le16();
    if ((len ^ 0xFFFF) != nlen) throw GzipError("corrupt stored block length");
    for (uint32_t k = 0; k < len; ++k) {
      room();
      put(static_cast<unsigned char>(byte()));
    }
  }

  void dynamic() {
    const int hlit = static_cast<int>(take(5)) + 257, hdist = static_cast<int>(take(5)) + 1;
    const int hclen = static_cast<int>(take(4)) + 4;
    if (hlit > 286 || hdist > 30) throw GzipError("invalid dynamic Huffman header");
    uint8_t cl[19] = {};
    for (int k = 0; k < hclen; ++k) cl[kClOrder[k]] = static_cast<uint8_t>(take(3));
    Huffman clh;
    clh.build(cl, 19);
    uint8_t len[286 + 30] = {};
    for (int k = 0; k < hlit + hdist;) {
      const int sym = decode(clh);
      if (sym < 16) {
        len[k++] = static_cast<uint8_t>(sym);
        continue;
      }
      int rep;
      uint8_t v = 0;
      if (sym == 16) {
        if (k == 0) throw GzipError("invalid code length repeat");
        v = len[k - 1];
        rep = 3 + static_cast<int>(take(2));
      } else if (sym == 17) {
        rep = 3 + static_cast<int>(take(3));
      } else {
        rep = 11 + static_cast<int>(take(7));
      }
      if (k + rep > hlit + hdist) throw GzipError("too many code lengths");
      while (rep-- > 0) len[k++] = v;
    }
    if (len[256] == 0) throw GzipError("missing end-of-block code");
    lit_.build(len, hlit);
    dist_.build(len + hlit, hdist);
    codes(lit_, dist_);
  }

  void codes(const Huffman& lit, const Huffman& dist) {
    for (;;) {
      room();
      const int sym = decode(lit);
      if (sym < 256) {
        put(static_cast<unsigned char>(sym));
        continue;
      }
      if (sym == 256) return;
      const int li = sym - 257;
      if (li >= 29) throw GzipError("invalid length code");
      const size_t length = kLenBase[li] + take(kLenExtra[li]);
      const int ds = decode(dist);
      if (ds >= 30) throw GzipError("invalid distance code");
      const size_t d = kDistBase[ds] + take(kDistExtra[ds]);
      if (d > pos_) throw GzipError("distance too far back");
      // Byte by byte: the source may overlap the destination (d < length repeats a pattern).
      unsigned char* o = out_.get();
      size_t src = pos_ - d;
      for (size_t k = 0; k < length; ++k) o[pos_ + k] = o[src + k];
      pos_ += length;
      member_ += length;
      total_ += length;
    }
  }

  const GzipSource& source_;
  const GzipSink& sink_;
  std::vector<unsigned char> in_;
  size_t in_pos_ = 0, in_len_ = 0;
  bool in_eof_ = false;
  uint64_t bitbuf_ = 0;
  int bitcnt_ = 0;
  std::unique_ptr<unsigned char[]> out_;
  size_t pos_ = 0, flushed_ = 0;
  uint32_t crc_ = 0;
  uint64_t member_ = 0, total_ = 0;
  Huffman fixed_lit_, fixed_dist_, lit_, dist_;
};

}  // namespace

uint32_t crc32_update(uint32_t crc, const unsigned char* data, size_t n) {
  const auto& t = crc_table();
  uint32_t c = crc ^ 0xFFFFFFFFu;
  for (size_t i = 0; i < n; ++i) c = t[(c ^ data[i]) & 0xFF] ^ (c >> 8);
  return c ^ 0xFFFFFFFFu;
}

bool is_gzip(const void* data, size_t n) {
  const auto* p = static_cast<const unsigned char*>(data);
  return n >= 2 && p[0] == 0x1F && p[1] == 0x8B;
}

uint64_t gunzip(const GzipSource& source, const GzipSink& sink) {
  Inflater inf(source, sink);
  return inf.run();
}

namespace {

std::atomic<uint64_t> g_max_output{kDefaultGzipMaxOutput};

// Appends to out, throwing once the total would exceed max_output.
void append_capped(std::string& out, const char* data, size_t n, uint64_t max_output) {
  if (n > max_output || out.size() > max_output - n) {
    throw GzipError("decompressed data larger than the limit of " + std::to_string(max_output) +
                    " bytes (corrupt or malicious gzip data? see set_gzip_max_output)");
  }
  out.append(data, n);
}

}  // namespace

void set_gzip_max_output(uint64_t bytes) { g_max_output.store(bytes == 0 ? kDefaultGzipMaxOutput : bytes); }

uint64_t gzip_max_output() { return g_max_output.load(); }

std::string gunzip(const std::string& compressed, uint64_t max_output) {
  size_t pos = 0;
  std::string out;
  const GzipSource src = [&](unsigned char* buf, size_t cap) {
    const size_t n = std::min(cap, compressed.size() - pos);
    std::memcpy(buf, compressed.data() + pos, n);
    pos += n;
    return n;
  };
  const GzipSink sink = [&](const char* data, size_t n) { append_capped(out, data, n, max_output); };
  gunzip(src, sink);
  return out;
}

std::string read_maybe_gzipped_file(const std::string& path, uint64_t max_output) {
  std::FILE* f = std::fopen(path.c_str(), "rb");
  if (!f) throw std::runtime_error(path + ": cannot open file");
  std::unique_ptr<std::FILE, int (*)(std::FILE*)> guard(f, &std::fclose);
  std::error_code ec;
  const auto size = std::filesystem::file_size(path, ec);
  unsigned char head[2];
  const size_t nh = std::fread(head, 1, 2, f);
  std::string data;
  if (is_gzip(head, nh)) {
    if (!ec) data.reserve(static_cast<size_t>(std::min<uintmax_t>({4 * size, uintmax_t{1} << 28, max_output})));
    bool head_sent = false;
    const GzipSource src = [&](unsigned char* buf, size_t cap) -> size_t {
      if (!head_sent) {
        head_sent = true;
        std::memcpy(buf, head, 2);
        return 2 + std::fread(buf + 2, 1, cap - 2, f);
      }
      const size_t got = std::fread(buf, 1, cap, f);
      if (got == 0 && std::ferror(f)) throw std::runtime_error(path + ": read error");
      return got;
    };
    const GzipSink sink = [&](const char* d, size_t n) { append_capped(data, d, n, max_output); };
    try {
      gunzip(src, sink);
    } catch (const GzipError& e) {
      throw GzipError(path + ": " + e.what());
    }
    return data;
  }
  if (!ec) data.reserve(static_cast<size_t>(size));
  data.append(reinterpret_cast<const char*>(head), nh);
  char buf[1 << 16];
  size_t got;
  while ((got = std::fread(buf, 1, sizeof(buf), f)) > 0) data.append(buf, got);
  if (std::ferror(f)) throw std::runtime_error(path + ": read error");
  return data;
}

}  // namespace milp
