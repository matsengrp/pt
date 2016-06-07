#include "partition.hpp"
#include "pll-utils.hpp"

/// @file partition.cpp
/// @brief Implementation of the Partition class.
/// Much of this was directly copied from libpll examples, so parts (c) Thomas Flouri.


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
  if(!tree_) {
    std::cout << "Parsing failure: tree must be an unrooted binary tree.\n";
    exit(EXIT_FAILURE);
  };

  if(!TreeHealthy(tree_))
    fatal("Missing branch lengths in tree.\n");

  char ** headers = NULL;
  char ** seqdata = NULL;
  unsigned int sites_count =
    ParseFasta(fasta_path, tip_nodes_count(), &headers, &seqdata);

  partition_ = pll_partition_create(
    tip_nodes_count(),      // Number of tip sequences we want to have.
    inner_nodes_count(),    // Number of CLV buffers to be allocated for inner nodes.
    STATES,                 // Number of states that our data have.
    sites_count,            // Number of sites in the alignment.
    1,                      // Number of different substitution models (or eigen decomposition)
                            //     to use concurrently (i.e. 4 for LG4).
    branch_count(),         // Number of probability matrices to be allocated.
    RATE_CATS,              // Number of rate categories we will use.
    inner_nodes_count(),    // How many scale buffers to use.
    ARCH_FLAGS);            // List of flags for hardware acceleration.

  SetModelParameters(partition_);

  EquipPartitionWithData(
    partition_, tree_, tip_nodes_count(), headers, seqdata);

  for(unsigned int i = 0; i < tip_nodes_count(); ++i) {
    free(seqdata[i]);
    free(headers[i]);
  }
  free(seqdata);
  free(headers);

  params_indices_ = (unsigned int *)malloc(RATE_CATS * sizeof(unsigned int));
  for(unsigned int i = 0; i < RATE_CATS; ++i) {
    params_indices_[i] = 0;
  }
  travbuffer_ = (pll_utree_t **)malloc(nodes_count() * sizeof(pll_utree_t *));
  branch_lengths_ = (double *)malloc(branch_count() * sizeof(double));
  matrix_indices_ = (unsigned int *)malloc(branch_count() * sizeof(unsigned int));
  operations_ = (pll_operation_t *)malloc(inner_nodes_count() *
                                                sizeof(pll_operation_t));

}


/// @brief Destructor for Partition.
Partition::~Partition() {
  free(params_indices_);
  free(travbuffer_);
  free(branch_lengths_);
  free(matrix_indices_);
  free(operations_);
  pll_partition_destroy(partition_);
  pll_utree_destroy(tree_);
}


/// @brief Returns a the tree as a Newick string.
/// @return Newick string.
std::string Partition::ToNewick() {
  return std::string(pll_utree_export_newick(tree_));
}


/// @brief Perform a full tree traversal and return the likelihood.
/// Eventually, we will want to break this up for efficiency gains
/// so that we aren't always doing a full traversal every time we
/// want to calculate the likelihood, but we'll do that carefully!
/// @return Log likelihood.
double Partition::FullTraversalLogLikelihood() {

  /* Perform a full postorder traversal of the unrooted tree. */
  unsigned int traversal_size;
  unsigned int matrix_count, ops_count;
  if(!pll_utree_traverse(
      tree_,
      cb_full_traversal,
      travbuffer_,
      &traversal_size)) {
    fatal("Function pll_utree_traverse() requires inner nodes as parameters");
  }

  /* Given the computed traversal descriptor, generate the operations
     structure, and the corresponding probability matrix indices that
     may need recomputing. */
  pll_utree_create_operations(
    travbuffer_,
    traversal_size,
    branch_lengths_,
    matrix_indices_,
    operations_,
    &matrix_count,
    &ops_count);

  /* Update matrix_count probability matrices for model with index 0. The i-th
     matrix (i ranges from 0 to matrix_count - 1) is generated using branch
     length branch_lengths[i] and can be referred to with index
     matrix_indices[i]. */
  pll_update_prob_matrices(
    partition_,
    params_indices_,
    matrix_indices_,
    branch_lengths_,
    matrix_count);

  /* Use the operations array to compute all ops_count inner CLVs. Operations
     will be carried out sequentially starting from operation 0 towrds ops_count-1. */
  pll_update_partials(partition_, operations_, ops_count);

  double logl = pll_compute_edge_loglikelihood(
    partition_,
    tree_->clv_index,
    tree_->scaler_index,
    tree_->back->clv_index,
    tree_->back->scaler_index,
    tree_->pmatrix_index,
    params_indices_,
    NULL); // Can supply a persite_lnl parameter.

  return logl;
}


}
