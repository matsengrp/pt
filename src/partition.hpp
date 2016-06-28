#ifndef PT_PARTITION_
#define PT_PARTITION_

#include <iostream>
#include <memory>
#include <pthread.h>
#include <string>

#include "libcuckoo/src/cuckoohash_map.hh"
#include "pll-utils.hpp"

/// @file partition.hpp
/// @brief Headers for the Partition class.

namespace pt {

/// @brief The representation of a tree, alignment, and all associated data.
///
class Partition {
  typedef cuckoohash_map<std::string, double> InnerTable;

private:
  unsigned int tip_nodes_count_;
  pll_utree_t *tree_;
  // Stores probability matrices, scalers, etc.
  pll_partition_t *partition_;
  unsigned int *matrix_indices_;
  double *branch_lengths_;
  pll_operation_t *operations_;
  // buffer for storing pointers to nodes of the tree in postorder traversal.
  pll_utree_t **travbuffer_;
  unsigned int *params_indices_;
  double *sumtable_;

public:
  Partition(std::string newick_path, std::string fasta_path,
            std::string RAxML_info_path);
  virtual ~Partition();

  unsigned int tip_nodes_count() { return tip_nodes_count_; };
  unsigned int inner_nodes_count() { return (tip_nodes_count() - 2); };
  unsigned int nodes_count() {
    return (tip_nodes_count() + inner_nodes_count());
  };
  unsigned int branch_count() { return (nodes_count() - 1); };

  std::string ToNewick(pll_utree_t *tree);
  void FullTraversalUpdate(pll_utree_t *tree);
  double LogLikelihood(pll_utree_t *tree);
  double FullTraversalLogLikelihood(pll_utree_t *tree);
  double OptimizeCurrentBranch(pll_utree_t *tree);
  void TreeBranchLengthsAux(pll_utree_t *tree);
  void TreeBranchLengths(pll_utree_t *tree);
  void FullBranchOpt(pll_utree_t *tree);
  pll_utree_t *SetNewickRoot(pll_utree_t *tree);
  bool SetLabelRoot(std::string label, pll_utree_t *tree, pll_utree_t **root);
  void RecursiveOrderedNewick(pll_utree_t *tree);
  std::string RootNewickRecursive(pll_utree_t *tree);
  pll_utree_t *ToOrderedNewick(pll_utree_t *tree);
  pll_utree_t *NNIUpdate(pll_utree_t *tree, int move_type);
  void MakeTables();
  char *utree_short_newick(pll_utree_t *root);
  static char *newick_utree_recurse(pll_utree_t *root);
  void NNITraverse(pll_utree_t *tree, double lambda);
  void NNIComputeEdge(pll_utree_t *tree, double lambda);
  InnerTable good_;
  InnerTable bad_;
};
}

#endif // PT_PARTITION_
