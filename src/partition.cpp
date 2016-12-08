#include "partition.hpp"

#include <memory>

#include "ordered-tree.hpp"
#include "pll-utils.hpp"

/// @file partition.cpp
/// @brief Implementation of the Partition class.
/// Much of this was directly copied from libpll examples, so parts (c) Thomas
/// Flouri.

namespace pt {

/// @brief Create a Partition object from a tree and an alignment.
///
/// This factory function (and private constructor) ensures that
/// Partition objects are managed via shared pointers. See the
/// comments in pt::Partition::NNIComputeEdge for details on why.
///
/// @param[in] newick_path
/// Path to a file with a Newick-format tree string.
/// @param[in] fasta_path
/// Path to a FASTA-formatted alignment.
/// @param[in] RAxML_info_path
/// Path to a RAxML info file.
std::shared_ptr<Partition> Partition::Create(const std::string &newick_path,
                                             const std::string &fasta_path,
                                             const std::string &RAxML_info_path)
{
  return std::shared_ptr<Partition>(new Partition(newick_path, fasta_path,
                                                  RAxML_info_path));
}

/// @brief Create a Partition object from an existing Partition and a tree.
///
/// This factory function (and private constructor) ensures that
/// Partition objects are managed via shared pointers. See the
/// comments in pt::Partition::NNIComputeEdge for details on why.
///
/// @param[in] obj
/// Pointer to partition object to copy.
/// @param[in] tree
/// Initial topology to use. (For full copy use obj->tree_)
std::shared_ptr<Partition> Partition::Create(const Partition &obj,
                                             pll_utree_t *tree)
{
  return std::shared_ptr<Partition>(new Partition(obj, tree));
}

/// @brief Private constructor for Partition from a tree and an alignment.
/// @param[in] newick_path
/// Path to a file with a Newick-format tree string.
/// @param[in] fasta_path
/// Path to a FASTA-formatted alignment.
/// @param[in] RAxML_info_path
/// Path to a RAxML info file.
Partition::Partition(std::string newick_path, std::string fasta_path,
                     std::string RAxML_info_path) {
  // Parse the unrooted binary tree in newick format, and store the number of
  // tip nodes in tip_nodes_count.
  tree_ = pll_utree_parse_newick(&newick_path[0], &tip_nodes_count_);
  if (!tree_) {
    std::cout << "Parsing failure: " << newick_path
              << " must be an unrooted binary tree.\n";
    exit(EXIT_FAILURE);
  };
  if (!TreeHealthy(tree_))
    fatal("Missing branch lengths in tree.\n");

  // Load in sequences and RAxML info.
  char **headers = NULL;
  char **seqdata = NULL;
  sites_count_ = ParseFasta(fasta_path, tip_nodes_count(), &headers, &seqdata);

  partition_ = CreatePartition();

  SetModelParameters(partition_, RAxML_info_path);

  EquipPartitionWithData(partition_, tree_, tip_nodes_count(), headers,
                         seqdata);

  for (unsigned int i = 0; i < tip_nodes_count(); ++i) {
    free(seqdata[i]);
    free(headers[i]);
  }
  free(seqdata);
  free(headers);

  // Allocate lots of memory for various operations.
  params_indices_ = (unsigned int *)malloc(RATE_CATS * sizeof(unsigned int));
  for (unsigned int i = 0; i < RATE_CATS; ++i) {
    params_indices_[i] = 0;
  }
  travbuffer_ = (pll_utree_t **)malloc(nodes_count() * sizeof(pll_utree_t *));
  branch_lengths_ = (double *)malloc(branch_count() * sizeof(double));
  matrix_indices_ =
      (unsigned int *)malloc(branch_count() * sizeof(unsigned int));
  operations_ =
      (pll_operation_t *)malloc(inner_nodes_count() * sizeof(pll_operation_t));
  sumtable_ = (double *)pll_aligned_alloc(
      partition_->sites * partition_->rate_cats * partition_->states_padded *
          sizeof(double),
      ALIGNMENT);
}

/// @brief Private copy constructor for Partition.
/// @param[in] obj
/// Pointer to partition object to copy.
/// @param[in] tree
/// Initial tree to use. (For full copy use obj->tree_)
Partition::Partition(const Partition &obj, pll_utree_t *tree) {
  tip_nodes_count_ = obj.tip_nodes_count_;
  tree_ = pll_utree_clone(tree);
  pll_utree_every(tree_, cb_copy_clv_traversal);
  sites_count_ = obj.sites_count_;

  partition_ = CreatePartition();

  for (unsigned int i = 0; i < partition_->tips + partition_->clv_buffers;
       ++i) {
    for (unsigned int j = 0;
         j < partition_->sites * partition_->states * partition_->rate_cats;
         j++)
      partition_->clv[i][j] = obj.partition_->clv[i][j];
  }

  for (unsigned int i = 0; i < partition_->prob_matrices; ++i) {
    for (unsigned int j = 0;
         j < partition_->states * partition_->states * partition_->rate_cats;
         j++)
      partition_->pmatrix[i][j] = obj.partition_->pmatrix[i][j];
  }

  for (unsigned int j = 0; j < partition_->rate_cats; j++)
    partition_->rates[j] = obj.partition_->rates[j];

  for (unsigned int i = 0;
       i < partition_->states * (partition_->states - 1) / 2; i++) {
    partition_->subst_params[0][i] = obj.partition_->subst_params[0][i];
  }

  for (unsigned int i = 0; i < partition_->rate_matrices; ++i) {
    for (unsigned int j = 0; j < partition_->states; j++)
      partition_->frequencies[i][j] = obj.partition_->frequencies[i][j];
  }

  for (unsigned int i = 0; i < partition_->rate_matrices; ++i) {
    for (unsigned int j = 0; j < partition_->states * partition_->states; j++)
      partition_->eigenvecs[i][j] = obj.partition_->eigenvecs[i][j];
  }

  for (unsigned int i = 0; i < partition_->rate_matrices; ++i) {
    for (unsigned int j = 0; j < partition_->states * partition_->states; j++)
      partition_->inv_eigenvecs[i][j] = obj.partition_->inv_eigenvecs[i][j];
  }

  for (unsigned int i = 0; i < partition_->rate_matrices; ++i) {
    for (unsigned int j = 0; j < partition_->states; j++)
      partition_->eigenvals[i][j] = obj.partition_->eigenvals[i][j];
  }

  partition_->maxstates = obj.partition_->maxstates;

  params_indices_ = (unsigned int *)malloc(RATE_CATS * sizeof(unsigned int));
  for (unsigned int i = 0; i < RATE_CATS; i++) {
    params_indices_[i] = obj.params_indices_[i];
  }

  travbuffer_ = (pll_utree_t **)malloc(nodes_count() * sizeof(pll_utree_t *));
  memcpy(travbuffer_, obj.travbuffer_, nodes_count() * sizeof(pll_utree_t *));

  branch_lengths_ = (double *)malloc(branch_count() * sizeof(double));
  memcpy(branch_lengths_, obj.branch_lengths_, branch_count() * sizeof(double));

  matrix_indices_ =
      (unsigned int *)malloc(branch_count() * sizeof(unsigned int));
  memcpy(matrix_indices_, obj.matrix_indices_,
         branch_count() * sizeof(unsigned int));

  operations_ =
      (pll_operation_t *)malloc(inner_nodes_count() * sizeof(pll_operation_t));
  memcpy(operations_, obj.operations_,
         inner_nodes_count() * sizeof(pll_operation_t));

  sumtable_ = (double *)pll_aligned_alloc(
      partition_->sites * partition_->rate_cats * partition_->states_padded *
          sizeof(double),
      ALIGNMENT);
  memcpy(sumtable_, obj.sumtable_, partition_->sites * partition_->rate_cats *
                                       partition_->states_padded *
                                       sizeof(double));
}

/// @brief Destructor for Partition.
Partition::~Partition() {
  free(params_indices_);
  free(travbuffer_);
  free(branch_lengths_);
  free(matrix_indices_);
  free(operations_);
  pll_aligned_free(sumtable_);
  pll_partition_destroy(partition_);
  pll_utree_every(tree_, cb_erase_data);
  pll_utree_destroy(tree_);
}

/// @brief Create a libpll partition data structure using the data in the
/// Partition object.
/// @return pll_partition_t data structure.
pll_partition_t *Partition::CreatePartition() {
  pll_partition_t *partition = pll_partition_create(
      tip_nodes_count(),   // Number of tip sequences we want to have.
      inner_nodes_count(), // Number of CLV buffers to be allocated for inner
                           // nodes.
      STATES,              // Number of states that our data have.
      sites_count_,        // Number of sites in the alignment.
      1, // Number of different substitution models (or eigen decomposition)
         //     to use concurrently (i.e. 4 for LG4).
      branch_count(),      // Number of probability matrices to be allocated.
      RATE_CATS,           // Number of rate categories we will use.
      inner_nodes_count(), // How many scale buffers to use.
      ARCH_FLAGS);         // List of flags for hardware acceleration.
  return partition;
}

/// @brief Perform a tree traversal and update CLV's, etc.
/// @param[in] tree
/// Tree.
/// @param[in] is_full
/// Which type of traversal update to perform. (1 = full) (0 = partial)
/// @return Number of nodes traversed.
unsigned int Partition::TraversalUpdate(pll_utree_t *tree, bool is_full) {
  unsigned int traversal_size;
  unsigned int matrix_count, ops_count;
  if (is_full) {
    if (!pll_utree_traverse(tree, cb_full_traversal, travbuffer_,
                            &traversal_size)) {
      fatal("Function pll_utree_traverse() requires inner nodes as parameters");
    }
  } else {
    if (!pll_utree_traverse(tree, cb_partial_traversal, travbuffer_,
                            &traversal_size)) {
      fatal("Function pll_utree_traverse() requires inner nodes as parameters");
    }
  }

  // Given the computed traversal descriptor, generate the operations
  // structure, and the corresponding probability matrix indices that
  // may need recomputing.
  pll_utree_create_operations(travbuffer_, traversal_size, branch_lengths_,
                              matrix_indices_, operations_, &matrix_count,
                              &ops_count);

  // Update matrix_count probability matrices for model with index 0. The i-th
  // matrix (i ranges from 0 to matrix_count - 1) is generated using branch
  // length branch_lengths[i] and can be referred to with index
  // matrix_indices[i].
  pll_update_prob_matrices(partition_, params_indices_, matrix_indices_,
                           branch_lengths_, matrix_count);

  // Use the operations array to compute all ops_count inner CLVs. Operations
  // will be carried out sequentially starting from operation 0 towrds
  // ops_count-1.
  pll_update_partials(partition_, operations_, ops_count);

  return traversal_size;
}

/// @brief Just calculate log likelihood.
/// NOTE: This will return nonsense if you haven't updated CLV's, etc.
/// @return Log likelihood.
double Partition::LogLikelihood(pll_utree_t *tree) {
  double logl = pll_compute_edge_loglikelihood(
      partition_, tree->clv_index, tree->scaler_index, tree->back->clv_index,
      tree->back->scaler_index, tree->pmatrix_index, params_indices_,
      nullptr); // Can supply a persite_lnl parameter as last argument.

  return logl;
}

/// @brief Make a traversal at a node and return the log likelihood.
/// @return Log likelihood.
double Partition::FullTraversalLogLikelihood(pll_utree_t *tree) {
  TraversalUpdate(tree, false);
  return LogLikelihood(tree);
}

/// @brief Optimize the current branch length, storing the new branch
/// length in the tree structure.
/// NOTE: This assumes we've just run `TraversalUpdate` or some such.
/// @param[in,out] tree
/// Parent node of branch to optimize.
/// @return Optimized branch length.
double Partition::OptimizeCurrentBranch(pll_utree_t *tree) {
  pll_utree_t *parent = tree;
  pll_utree_t *child = tree->back;

  // Compute the sumtable for the particular branch once before proceeding with
  // the optimization.
  pll_update_sumtable(partition_, parent->clv_index, child->clv_index,
                      params_indices_, sumtable_);

  double len = tree->length;
  bool maybe_decreasing = false;

  for (unsigned int i = 0; i < MAX_ITER; ++i) {
    double d1; // First derivative.
    double d2; // Second derivative.

    pll_compute_likelihood_derivatives(
        partition_, parent->scaler_index, child->scaler_index, len,
        params_indices_, sumtable_, &d1, &d2);

    // printf("Branch length: %f log-L: %f Derivative: %f D2: %f\n", len,
    // opt_logl, d1,d2);
    // If derivative is approximately zero then we've found the maximum.
    if (fabs(d1) < EPSILON)
      break;

    // Newton's method for finding the optimum of a function. The iteration to
    // reach the optimum is

    // x_{i+1} = x_i - f'(x_i) / f''(x_i)

    // where x_i is the current branch, f'(x_i) is the first derivative and
    // f''(x_i) is the second derivative of the likelihood function.
    if (d2 < 0)
      len += d1 / d2;
    else
      len -= d1 / d2;

    // If the next branch length to evaluate goes negative, we instead
    // set it to a small positive value for the next iteration. If
    // this has happened before, we stop early, as the curve is
    // probably decreasing.
    if (len < 0) {
      len = EPSILON;

      if (maybe_decreasing) {
        break;
      }

      maybe_decreasing = true;
    }
  }

  // Update current branch lengths.
  parent->length = len;
  child->length = len;

  // Update this branch's probability matrix now that the branch
  // length has changed. No CLVs need to be invalidated.
  pll_update_prob_matrices(partition_, params_indices_, &(parent->pmatrix_index),
                           &(parent->length), 1);

  return len;
}

/// @brief Aux function for optimizing branch lengths across the whole tree.
/// @param[in] tree
/// Child node from which to optimize the branch length and continue recursion.
void Partition::TreeBranchLengthsAux(pll_utree_t *tree) {
  if (!tree->next) {
    TraversalUpdate(tree->back, false);
    OptimizeCurrentBranch(tree->back);
  } else {
    TraversalUpdate(tree, false);
    OptimizeCurrentBranch(tree);
    TreeBranchLengthsAux(tree->next->back);
    TreeBranchLengthsAux(tree->next->next->back);
  }
}

///@brief Perform a postorder tree traversal which optimizes the branch length
/// at every edge.
/// NOTE: This optimizes the initial branch twice per traversal.
/// @param[in] tree
/// Tree.
void Partition::TreeBranchLengths(pll_utree_t *tree) {
  if (!tree->next) {
    fatal("Function TreeBranchLengths requires an inner node as parameter");
  }
  TreeBranchLengthsAux(tree);
  TreeBranchLengthsAux(tree->back);
}

///@brief Repeat branch length optimization until log likelihood changes very
/// little or MAX_ITER is reached.
void Partition::FullBranchOpt(pll_utree_t *tree) {
  double loglike_prev = FullTraversalLogLikelihood(tree);

  TreeBranchLengths(tree);

  double loglike = FullTraversalLogLikelihood(tree);

  unsigned int i = 0;
  while (fabs(loglike_prev - loglike) > EPSILON && i < MAX_ITER) {
    TreeBranchLengths(tree);

    loglike_prev = loglike;
    loglike = FullTraversalLogLikelihood(tree);
    i++;
  }
}

///@brief Perform an NNI move at the current edge and
/// order the tree.
///@param[in] tree
/// Node at edge on which to perform NNI.
///@param[in] move_type
/// Type of NNI move to perform.
///@return Ordered new topology.
pll_utree_t *Partition::NNIUpdate(pll_utree_t *tree, int move_type) {
  // Orient CLV's
  TraversalUpdate(tree, false);
  pll_utree_nni(tree, move_type, nullptr);
  // Recalculate CLV's after NNI
  node_info_t *node_info;
  node_info = (node_info_t *)tree->data;
  if (node_info) {
    node_info->clv_valid = 0;
  }
  if (tree->back->next) {
    node_info = (node_info_t *)tree->back->data;
    if (node_info) {
      node_info->clv_valid = 0;
    }
  }
  TraversalUpdate(tree, false);
  // Reorder the tree
  tree = ToOrderedNewick(tree);
  return tree;
}

/// @brief Add a MakeTables job to a thread pool work queue.
///
/// The caller guarantees that the arguments passed by pointer or
/// reference will continue to exist until all threads are completed.
/// See Partition::MakeTables for argument descriptions.
///
/// TODO: At time of writing, calls to MakeTables are made using the
/// tree_ member variable. As such the need to pass in a tree here
/// could likely be eliminated.
void Partition::QueueMakeTables(double cutoff, double logl, pll_utree_t *tree,
                                TreeTable &good, TreeTable &all,
                                ctpl::thread_pool &pool)
{
  auto this_shared = shared_from_this();
  pool.push([this_shared, cutoff, logl, tree, &good, &all, &pool](int) {
      this_shared->MakeTables(cutoff, logl, tree, good, all, pool);
    });
}

///@brief Perform every possible NNI move from the current state, and sort them
/// into good and all tables.
/// NOTE: This function assumes that the current topology is the ML tree.
/// @param[in] cutoff
/// Scaler value by which to multiply ML tree Log-L, the result is the cutoff
/// between good and bad trees.
/// @param[in] logl
/// The log likelihood of the ML tree.
/// @param[in] tree
/// Internal node of topology on which to try NNI moves.
/// @param[in] good
/// Hash table for good trees.
/// @param[in] all
/// Hash table for all trees.
/// @param[in] pool
/// Thread pool to which to push.
void Partition::MakeTables(double cutoff, double logl, pll_utree_t *tree,
                           TreeTable &good, TreeTable &all,
                           ctpl::thread_pool &pool) {
  // Order the tree, check if it is in the good table, and add if it is not.
  tree = ToOrderedNewick(tree);
  /* if (!good.contains(ToNewick(tree))) {
     good.insert(ToNewick(tree), logl);
     all.insert(ToNewick(tree), 0);
   }*/
  // Traverse the tree, performing both possible NNI moves, and sorting into
  // tables at each internal edge.
  NNITraverse(tree, logl, cutoff, good, all, pool);
  NNITraverse(tree->back, logl, cutoff, good, all, pool);
}

/// @brief Print out hashtables.
/// @param[in] print_all
/// Prints the "all trees" table if true.
/// @param[in] good
/// Table of good trees.
/// @param[in] all
/// Table of all trees.
void Partition::PrintTables(bool print_all, TreeTable &good, TreeTable &all) {
  // Print Tables.
  std::cout << "Good: " << std::endl;

  auto lt = good.lock_table();
  for (auto &item : lt)
    std::cout << "Log Likelihood for " << item.first << " : " << item.second
              << std::endl;

  if (print_all) {
    std::cout << "All: " << std::endl;

    auto lt1 = all.lock_table();
    for (auto &item : lt1)
      std::cout << item.first << std::endl;
  }
  printf("Trees Explored: %zu \n", all.size());
  printf("Good Trees Found: %zu \n", good.size());
}

/// @brief Perform both NNI moves at an edge and compare their log-likelihoods
/// to the ML, sort accordingly.
/// @param[in] tree
/// The current edge.
/// @param[in] lambda
/// The likelihood of ML tree.
/// @param[in] cutoff
/// The scaler cutoff for tree acceptance (c * lambda)
/// @param[in] good
/// Hash table for good trees.
/// @param[in] all
/// Hash table for all trees.
/// @param[in] pool
/// Thread pool to which to push jobs.
void Partition::NNIComputeEdge(pll_utree_t *tree, int move_type, double lambda,
                               double cutoff, TreeTable &good, TreeTable &all,
                               ctpl::thread_pool &pool) {
  // Create a clone of the original tree to perform NNI and reordering on.
  pll_utree_t *clone = pll_utree_clone(tree);
  // Deep copy clv_valid values.
  pll_utree_every(clone, cb_copy_clv_traversal);
  // Perform NNI and reordering on edge.
  clone = NNIUpdate(clone, move_type);

  // Keep track of whether or not we should free the cloned tree at
  // the end of this function, or if it will be needed later.
  bool free_clone = true;

  std::string label = ToNewick(clone);
  if (all.insert(label, 0)) {
    FullBranchOpt(clone);
    /// @todo Can we save time by having FullBranchOpt return this?
    double lambda_1 = FullTraversalLogLikelihood(clone);
    // Compare new likelihood to ML, then decide which table to put in.
    if (lambda_1 > cutoff * lambda) {
      if (good.insert(label, lambda_1)) {
        // Don't free the cloned tree at the end of this function.
        // Instead it will be freed by the function pushed to the
        // thread pool.
        free_clone = false;

        // Create routine for good tree and have it MakeTables. Push
        // routine to pool.
        //
        // This Partition object should only be destroyed after all
        // the jobs it adds to the pool are completed. Here we create
        // a shared pointer to this Partition to pass to the lambda to
        // ensure that it still exists when the new Partition inside
        // the lambda is created. This is only safe if this object is
        // already owned by a shared pointer, which is why the
        // constructors are private and object creation is delegated
        // to factory functions.
        auto this_shared = shared_from_this();
        pool.push([this_shared, clone, cutoff, lambda, &good, &all, &pool](int id) {
          auto temp = pt::Partition::Create(*this_shared, clone);
          temp->MakeTables(cutoff, lambda, temp->tree_, good, all, pool);
          pll_utree_every(clone, cb_erase_data);
          pll_utree_destroy(clone);
        });

        std::cout << "work queue size: " << pool.queue_size() << "\n";
      }
    }
  }

  if (free_clone) {
    pll_utree_every(clone, cb_erase_data);
    pll_utree_destroy(clone);
  }
}

/// @brief Traverse the tree and perform NNI moves at each internal edge.
void Partition::NNITraverse(pll_utree_t *tree, double lambda, double cutoff,
                            TreeTable &good, TreeTable &all,
                            ctpl::thread_pool &pool) {
  if (!tree->next)
    return;
  if (!tree->back->next) {
    NNITraverse(tree->next, lambda, cutoff, good, all, pool);
  }
  NNIComputeEdge(tree, PLL_UTREE_MOVE_NNI_LEFT, lambda, cutoff, good, all, pool);
  NNIComputeEdge(tree, PLL_UTREE_MOVE_NNI_RIGHT, lambda, cutoff, good, all, pool);
  NNITraverse(tree->next->back, lambda, cutoff, good, all, pool);
  NNITraverse(tree->next->next->back, lambda, cutoff, good, all, pool);
}
}
