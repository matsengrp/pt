#include "partition.hpp"

/// @file partition.cpp
/// @brief Implementation of the Partition class.

namespace pt {

#define STATES    4
#define RATE_CATS 4


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
  if(!tree_) {
    std::cout << "Parsing failure: tree must be an unrooted binary tree.\n";
    exit(EXIT_FAILURE);
  };

  /*
  if(!TreeHealthy(tree_)) {
    std::cout << "Missing branch lengths in tree.\n";
    exit(EXIT_FAILURE);
  }
  */

}


/// @brief Destructor for Partition.
Partition::~Partition() {

  pll_utree_destroy(tree_);

}


/// @brief Returns a the tree as a Newick string.
/// @return Newick string.
std::string Partition::ToNewick() {

  return std::string(pll_utree_export_newick(tree_));

}

}
