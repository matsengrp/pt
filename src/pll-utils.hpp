#ifndef PT_PLL_UTILS_
#define PT_PLL_UTILS_

#include <libpll/pll.h>
#include <string>

/// @file pll-utils.hpp
/// @brief Utilities for interfacing with libpll.
/// Much of this was directly copied from libpll examples, so parts (c) Thomas
/// Flouri.


namespace pt {

const unsigned int STATES = 4;
const unsigned int RATE_CATS = 4;
const unsigned int ARCH_FLAGS = PLL_ATTRIB_ARCH_AVX;

void fatal(const char* format, ...);
int cb_full_traversal(pll_utree_t* node);

bool TreeHealthy(pll_utree_t* tree);
unsigned int ParseFasta(std::string path, unsigned int seq_count,
                        char*** headers_out, char*** seqdata_out);
void EquipPartitionWithData(pll_partition_t* partition, pll_utree_t* tree,
                            unsigned int tip_nodes_count, char** headers,
                            char** seqdata);
void SetModelParameters(pll_partition_t* partition);
}

#endif  // PT_PLL_UTILS_
