#pragma once
// Resident memory of this process, for enforcing --mem-limit on everything the
// solver allocates (estimates of individual structures miss some of it).

namespace milp {

// Memory charged to this process in bytes (resident set; macOS: physical footprint), or 0 where it
// cannot be measured.
long long resident_bytes();

// Asks the allocator to return freed memory it keeps cached to the OS (macOS: malloc pressure
// relief; glibc: malloc_trim). A no-op elsewhere.
void release_free_memory();

}  // namespace milp
