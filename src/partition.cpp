#include "partition.hpp"

/// @file partition.cpp
/// @brief Implementation of the Partition class.

namespace pt {


/// @brief Constructor for Partition from a tree and an alignment.
/// @param[in] newick_path
/// Path to a file with a Newick-format tree string.
/// @param[in] fasta_path
/// Path to a FASTA-formatted alignment.
Partition::Partition(
  std::string newick_path,
  std::string fasta_path) {

  // Parse the unrooted binary tree in newick format, and store the number of
  // tip nodes in tip_nodes_count.
  tree_ = pll_utree_parse_newick(&newick_path[0], &tip_nodes_count_);
  if (!tree_) {
    std::cout << "Parsing failure: tree must be an unrooted binary tree.";
    exit(EXIT_FAILURE);
  };


}

Partition::~Partition() {

  pll_utree_destroy(tree_);

}
}
