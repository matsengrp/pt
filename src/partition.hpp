#ifndef PT_PARTITION_
#define PT_PARTITION_

#include <iostream>
#include <memory>

#include <cuckoohash_map.hh>

#include "ctpl_stl.h"
#include "pll-utils.hpp"

/// @file partition.hpp
/// @brief Headers for the Partition class.

namespace pt {

typedef cuckoohash_map<std::string, double> TreeTable;

/// @brief The representation of a tree, alignment, and all associated data.
class Partition : public std::enable_shared_from_this<Partition> {
private:
  pll_partition_t *partition_;
  unsigned int sites_count_;
  unsigned int tip_nodes_count_;
  unsigned int *matrix_indices_;
  double *branch_lengths_;
  pll_operation_t *operations_;
  // Buffer for storing pointers to nodes of the tree in postorder traversal.
  pll_utree_t **travbuffer_;
  unsigned int *params_indices_;
  double *sumtable_;

 public:
  static std::shared_ptr<Partition> Create(const std::string &newick_path,
                                           const std::string &fasta_path,
                                           const std::string &RAxML_info_path);
  static std::shared_ptr<Partition> Create(const Partition &obj,
                                           pll_utree_t *tree);

  virtual ~Partition();

  pll_utree_t *tree_;
  pll_partition_t* GetPartition() {return partition_ ; }
  const pll_partition_t *GetPartition() const { return partition_ ; }

  unsigned int tip_nodes_count() { return tip_nodes_count_; };
  unsigned int inner_nodes_count() { return (tip_nodes_count() - 2); };
  unsigned int nodes_count() {
    return (tip_nodes_count() + inner_nodes_count());
  };
  unsigned int branch_count() { return (nodes_count() - 1); };

  unsigned int TraversalUpdate(pll_utree_t *tree, TraversalType type);
  double LogLikelihood(pll_utree_t *tree);
  double FullTraversalLogLikelihood(pll_utree_t *tree);
  double OptimizeCurrentBranch(pll_utree_t *tree);
  void FullBranchOpt(pll_utree_t *tree);
  void QueueMakeTables(double cutoff, double logl, pll_utree_t *tree,
                       TreeTable &good, TreeTable &all, ctpl::thread_pool &pool);
  void PrintTables(bool print_all, TreeTable &good, TreeTable &all);

 private:
  Partition(std::string newick_path, std::string fasta_path,
            std::string RAxML_info_path);
  Partition(const Partition &obj, pll_utree_t *tree);

  pll_partition_t *CreatePartition();
  void TreeBranchLengthsAux(pll_utree_t *tree);
  void TreeBranchLengths(pll_utree_t *tree);
  pll_utree_t *NNIUpdate(pll_utree_t *tree, int move_type);
  void NNITraverse(pll_utree_t *tree, double lambda, double cutoff,
                   TreeTable &good, TreeTable &all, ctpl::thread_pool &pool);
  void NNIComputeEdge(pll_utree_t *tree, int move_type, double lambda,
                      double cutoff, TreeTable &good, TreeTable &all,
                      ctpl::thread_pool &pool);
  void MakeTables(double cutoff, double logl, pll_utree_t *tree,
                  TreeTable &good, TreeTable &all, ctpl::thread_pool &pool);
};
}

#endif // PT_PARTITION_
