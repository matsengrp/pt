#ifndef PT_PLL_UTILS_
#define PT_PLL_UTILS_

#include <string>
#include <vector>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

/// @file pll-utils.hpp
/// @brief Utilities for interfacing with libpll.
/// Much of this was directly copied from libpll examples, so parts (c) Thomas
/// Flouri.

namespace pt {

typedef struct { int clv_valid; } node_info_t;

const unsigned int STATES = 4;
const unsigned int RATE_CATS = 4;
const unsigned int ALIGNMENT = PLL_ALIGNMENT_AVX;
const unsigned int ARCH_FLAGS = PLL_ATTRIB_ARCH_AVX;
// MAX_ITER: Max iterations when optimizing branch lengths.
const unsigned int MAX_ITER = 32;
const double EPSILON = 1e-5; // Threshold for detecting zero.

enum class TraversalType { FULL, PARTIAL };

int cb_full_traversal(pll_utree_t *node);
int cb_partial_traversal(pll_utree_t *node);
int cb_copy_clv_traversal(pll_utree_t *node);
int cb_branch_healthy(pll_utree_t *tree);
int cb_erase_data(pll_utree_t *tree);

pll_partition_t *pllext_partition_clone(pll_partition_t *rhs);

bool TreeHealthy(pll_utree_t *tree);
unsigned int ParseFasta(std::string path, unsigned int seq_count,
                        std::vector<std::string> &headers_out,
                        std::vector<std::string> &seqdata_out);
void EquipPartitionWithData(pll_partition_t *partition, pll_utree_t *tree,
                            unsigned int tip_nodes_count,
                            const std::vector<std::string>& headers,
                            const std::vector<std::string>& seqdata);
std::vector<std::string> ssplit(const std::string &s, char delim);
void SetModelParameters(pll_partition_t *partition, std::string path);
}

#endif // PT_PLL_UTILS_
