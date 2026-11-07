#pragma once
// Anonymized copies of a model (popso --write-anon), for users who cannot share a model as it is
// but can share one without its meaning: generic names and a shuffled row and column order, with
// every number unchanged, so the copy reproduces the solver's behaviour up to the order (some
// bugs depend on the order; AnonOptions::keep_order keeps it).
//
// The copy has columns x1..xn, rows c1..cm, SOS sets s1.., indicator constraints i1.. and the
// model name "anon". Columns, rows, SOS sets and indicator constraints are each permuted by a
// seeded shuffle (the same seed gives the same copy on every platform: its own Fisher-Yates with
// SplitMix64, not std::shuffle), and every part follows the permutations: the matrix, bounds,
// types, objective, quadratic objective (lower-triangle convention kept), SOS members and
// weights, indicator columns, coefficients and binaries. Within a column (matrix), a quadratic
// objective and an indicator row, entries are sorted by the new indices, so that their order
// does not reveal the original one. write_mps writes the copy with round-trip precision (up to
// the ranged rows MPS cannot express exactly, which it counts; see mps_reader.hpp).

#include <cstdint>
#include <string>
#include <vector>

#include "milp/core/model.hpp"

namespace milp {

inline constexpr uint64_t kAnonDefaultSeed = 1;

struct AnonOptions {
  uint64_t seed = kAnonDefaultSeed;
  bool keep_order = false;  // no shuffle: only the names change
};

// The permutations of an anonymization: entry k of the copy is entry col[k] (row[k], sos[k],
// indicator[k]) of the original.
struct AnonMap {
  std::vector<Index> col, row, sos, indicator;
};

Model anonymize(const Model& model, const AnonOptions& options = {}, AnonMap* map = nullptr);

// The private mapping back to the original names (popso --anon-map): one line
// "<new name> <original name>" per column, row, SOS set and indicator constraint of the copy, in
// that order (columns x1.., rows c1.., sets s1.., indicators i1..), after a '#' comment line.
// An unnamed original appears as C<j>, R<i>, S<k> or I<k> (its 0-based index); whitespace in a
// name becomes '_' (as write_mps writes it), so every line has exactly two fields.
std::string anon_map_string(const Model& original, const AnonMap& map);
void write_anon_map(const std::string& path, const Model& original, const AnonMap& map);

}  // namespace milp
