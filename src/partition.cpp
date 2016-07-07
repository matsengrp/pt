#include "partition.hpp"
#include "pll-utils.hpp"
#include <thread>

/// @file partition.cpp
/// @brief Implementation of the Partition class.
/// Much of this was directly copied from libpll examples, so parts (c) Thomas
/// Flouri.

namespace pt {
/// @brief Constructor for Partition from a tree and an alignment.
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

  char **headers = NULL;
  char **seqdata = NULL;
  unsigned int sites_count =
      ParseFasta(fasta_path, tip_nodes_count(), &headers, &seqdata);

  partition_ = pll_partition_create(
      tip_nodes_count(),   // Number of tip sequences we want to have.
      inner_nodes_count(), // Number of CLV buffers to be allocated for inner
                           // nodes.
      STATES,              // Number of states that our data have.
      sites_count,         // Number of sites in the alignment.
      1, // Number of different substitution models (or eigen decomposition)
         //     to use concurrently (i.e. 4 for LG4).
      branch_count(),      // Number of probability matrices to be allocated.
      RATE_CATS,           // Number of rate categories we will use.
      inner_nodes_count(), // How many scale buffers to use.
      ARCH_FLAGS);         // List of flags for hardware acceleration.

  SetModelParameters(partition_, RAxML_info_path);

  EquipPartitionWithData(partition_, tree_, tip_nodes_count(), headers,
                         seqdata);

  for (unsigned int i = 0; i < tip_nodes_count(); ++i) {
    free(seqdata[i]);
    free(headers[i]);
  }
  free(seqdata);
  free(headers);

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

/// @brief Destructor for Partition.
Partition::~Partition() {
  free(params_indices_);
  free(travbuffer_);
  free(branch_lengths_);
  free(matrix_indices_);
  free(operations_);
  pll_aligned_free(sumtable_);
  pll_partition_destroy(partition_);
  pll_utree_destroy(tree_);
}

/// @brief Returns a the tree as a Newick string.
/// @return Newick string.
std::string Partition::ToNewick(pll_utree_t *tree) {
  return std::string(utree_short_newick(tree));
}

/// @brief Recursive function to generate newick string without branch lengths.
char *Partition::newick_utree_recurse(pll_utree_t *root) {
  char *newick;
  int size_alloced;
  assert(root != NULL);
  if (!root->next)
    size_alloced = asprintf(&newick, "%s", root->label);
  else {
    char *subtree1 = newick_utree_recurse(root->next->back);
    if (subtree1 == NULL)
      return NULL;
    char *subtree2 = newick_utree_recurse(root->next->next->back);
    if (subtree2 == NULL) {
      free(subtree1);
      return NULL;
    }
    size_alloced = asprintf(&newick, "(%s,%s)%s", subtree1, subtree2,
                            root->label ? root->label : "");
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0) {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
    return NULL;
  }

  return newick;
}

/// @brief Prints newick string without branch lengths.
char *Partition::utree_short_newick(pll_utree_t *root) {
  char *newick;
  int size_alloced;
  if (!root)
    return NULL;

  if (!root->next)
    root = root->back;

  char *subtree1 = newick_utree_recurse(root->back);
  if (subtree1 == NULL)
    return NULL;
  char *subtree2 = newick_utree_recurse(root->next->back);
  if (subtree2 == NULL) {
    free(subtree1);
    return NULL;
  }
  char *subtree3 = newick_utree_recurse(root->next->next->back);
  if (subtree3 == NULL) {
    free(subtree1);
    free(subtree2);
    return NULL;
  }

  size_alloced = asprintf(&newick, "(%s,%s,%s)%s;", subtree1, subtree2,
                          subtree3, root->label ? root->label : "");
  free(subtree1);
  free(subtree2);
  free(subtree3);
  if (size_alloced < 0) {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
    return NULL;
  }

  return (newick);
}

/// @brief Finds node in tree's subtree with the label that is first
/// alphabetically.
/// @param[in] tree
/// Parent node of subtree to root.
/// @return smallest label
std::string Partition::FindRootNode(pll_utree_t *tree) {
  std::string minlabel;
  if (!tree->next)
    return tree->label;
  std::string rlabel1 = FindRootNode(tree->next->back);
  std::string rlabel2 = FindRootNode(tree->next->next->back);
  if (rlabel1.compare(rlabel2) < 0)
    minlabel = rlabel1;
  else {
    minlabel = rlabel2;
    pll_utree_t *temp = tree->next;
    tree->next = tree->next->next;
    tree->next->next = temp;
    tree->next->next->next = tree;
  }

  return minlabel;
}

/// @brief Finds a node with the given label and sets the current node to the
/// parent of it.
/// @param[in] label
/// Label of node to find
/// @param[in] tree
/// Node to search at
/// @param[in] root
/// Buffer for storing root node
bool Partition::SetLabelRoot(std::string label, pll_utree_t *tree,
                             pll_utree_t **root) {
  if (!tree->next) {
    if (tree->label == label) {
      root[0] = tree->back;
      return true;
    } else
      return false;
  } else {
    return (SetLabelRoot(label, tree->next->back, root) ||
            SetLabelRoot(label, tree->next->next->back, root));
  }
}

/// @brief Sets the label of the root for the entire tree.
/// @param[in] tree
/// An internal node.
/// @return Parent of the first node alphabetically.
pll_utree_t *Partition::SetNewickRoot(pll_utree_t *tree) {
  std::string minlabel;
  std::string rlabel1 = FindRootNode(tree);
  std::string rlabel2 = FindRootNode(tree->back);
  if (rlabel1.compare(rlabel2) < 0)
    minlabel = rlabel1;
  else
    minlabel = rlabel2;
  pll_utree_t **root;
  root = (pll_utree_t **)malloc(sizeof(pll_utree_t *));
  SetLabelRoot(minlabel, tree, root);
  SetLabelRoot(minlabel, tree->back, root);
  return root[0];
}

/// @brief Determines order for first ternary step, then runs recursive ordering
/// for the two non-root subtrees.
/// @return Completely ordered newick string.
pll_utree_t *Partition::ToOrderedNewick(pll_utree_t *tree) {
  tree = SetNewickRoot(tree);
  FindRootNode(tree);
  return tree;
}

/// @brief Perform either a full or fast traversal update.
/// @param[in] tree
/// Parent node from which to update CLV's etc.
/// @param[in] is_full
/// Which type of traversal update to perform. (1 = full)
void Partition::TraversalUpdate(pll_utree_t *tree, bool is_full) {
  // Perform a full postorder traversal of the unrooted tree.
  if (is_full)
    FullTraversalUpdate(tree);
  else
    FastUpdate(tree);
}

/// @brief Perform a full tree traversal and update CLV's, etc.
/// Eventually, we will want to break this up for efficiency gains
/// so that we aren't always doing a full traversal every time we
/// want to calculate the likelihood, but we'll do that carefully!
/// @param[in] tree
/// Parent node from which to update CLV's etc.
void Partition::FullTraversalUpdate(pll_utree_t *tree) {
  unsigned int traversal_size;
  unsigned int matrix_count, ops_count;
  if (!pll_utree_traverse(tree, cb_full_traversal, travbuffer_,
                          &traversal_size)) {
    fatal("Function pll_utree_traverse() requires inner nodes as parameters");
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
}

void Partition::FastUpdate(pll_utree_t *tree) {
  // Some other type of updating here.
}

/// @brief Just calculate log likelihood.
/// NOTE: This will return nonsense if you haven't updated CLV's, etc.
/// @return Log likelihood.
double Partition::LogLikelihood(pll_utree_t *tree) {
  double logl = pll_compute_edge_loglikelihood(
      partition_, tree->clv_index, tree->scaler_index, tree->back->clv_index,
      tree->back->scaler_index, tree->pmatrix_index, params_indices_,
      NULL); // Can supply a persite_lnl parameter as last argument.

  return logl;
}

/// @brief Make a traversal at a node and return the log likelihood.
/// @return Log likelihood.
double Partition::FullTraversalLogLikelihood(pll_utree_t *tree) {
  TraversalUpdate(tree, 1);
  return LogLikelihood(tree);
}

/// @brief Optimize the current branch length.
/// NOTE: This assumes we've just run `FullTraversalUpdate` or some such.
/// @param[in] tree
/// Parent node of branch to optimize.
/// @return Optimized branch length.
double Partition::OptimizeCurrentBranch(pll_utree_t *tree) {
  // Perform a full postorder traversal of the unrooted tree.
  pll_utree_t *parent = tree;
  pll_utree_t *child = tree->back;

  // Compute the sumtable for the particular branch once before proceeding with
  // the optimization.
  pll_update_sumtable(partition_, parent->clv_index, child->clv_index,
                      params_indices_, sumtable_);

  double len = tree->length;
  double d1; // First derivative.
  double d2; // Second derivative.
  for (unsigned int i = 0; i < MAX_ITER; ++i) {
    double opt_logl = pll_compute_likelihood_derivatives(
        partition_, parent->scaler_index, child->scaler_index, len,
        params_indices_, sumtable_, &d1, &d2);

    /// printf("Branch length: %f log-L: %f Derivative: %f D2: %f\n", len,
    /// opt_logl, d1,d2);
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

    // Branch optimization was returning negative values, set minimum branch
    // length to be small positive value

    if (len < 0) {
      len = EPSILON;
      double opt_logl = pll_compute_likelihood_derivatives(
          partition_, parent->scaler_index, child->scaler_index, len,
          params_indices_, sumtable_, &d1, &d2);
      /// printf("Branch length: %f log-L: %f Derivative: %f D2: %f\n", len,
      /// opt_logl, d1,d2);
      break;
    }
  }

  // Update current branch lengths.
  parent->length = len;
  child->length = len;

  return len;
}

// Potential callback function for branch length optimization.
/*int cb_branch_lengths(pll_utree_t *tree){
  if (!tree->next) {
    TraversalUpdate(tree->back, 1);
    OptimizeCurrentBranch(tree->back);
  } else {
    TraversalUpdate(tree, 1);
    OptimizeCurrentBranch(tree);
  }
  return 1;
}*/

/// @brief Aux function for optimizing branch lengths across the whole tree.
/// @param[in] tree
/// Child node from which to optimize the branch length and continue recursion.
void Partition::TreeBranchLengthsAux(pll_utree_t *tree) {
  if (!tree->next) {
    TraversalUpdate(tree->back, 1);
    OptimizeCurrentBranch(tree->back);
  } else {
    TraversalUpdate(tree, 1);
    OptimizeCurrentBranch(tree);
    TreeBranchLengthsAux(tree->next->back);
    TreeBranchLengthsAux(tree->next->next->back);
  }
}

///@brief Perform a postorder tree traversal which optimizes the branch length
/// at every edge.
/// NOTE: This optimizes the initial branch twice per traversal.
void Partition::TreeBranchLengths(pll_utree_t *tree) {
  if (!tree->next) {
    fatal("Function TreeBranchLengthsAux requires an inner node as parameter");
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

  int i = 0;
  while (fabs(loglike_prev - loglike) > EPSILON && i < MAX_ITER) {
    TreeBranchLengths(tree);

    loglike_prev = loglike;
    loglike = FullTraversalLogLikelihood(tree);
    i++;
  }
}

///@brief Perform an NNI move at the current edge, optimize branch lengths, and
/// order the tree.
///@param[in] tree
/// Node at edge on which to perform NNI.
///@param[in] move_type
/// Type of NNI move to perform.
///@return Ordered, optimized new topology.
pll_utree_t *Partition::NNIUpdate(pll_utree_t *tree, int move_type) {
  pll_utree_rb_t *rb;
  rb = (pll_utree_rb_t *)malloc(sizeof(pll_utree_t *) + sizeof(int));
  pll_utree_nni(tree, move_type, rb);
  tree = ToOrderedNewick(tree);
  return tree;
}

///@brief Perform every possible NNI move from the current state, and sort them
/// into good and bad tables.
/// NOTE: This function assumes that the current topology is "good" (i.e. it is
/// the ML tree).
void Partition::MakeTables(double cutoff, double logl, pll_utree_t *tree) {
  // Update and optimize the ML tree, store its logl for comparison, and add it
  // to the good table.
  pll_utree_t *clone = pll_utree_clone(tree);
  clone = ToOrderedNewick(clone);
  if (!good_.contains(ToNewick(clone))) {
    good_.insert(ToNewick(clone), logl);
  }
  // Traverse the tree, performing both possible NNI moves, and sorting into
  // tables at each internal edge.
  NNITraverse(tree, logl, cutoff);
  NNITraverse(tree->back, logl, cutoff);
}
void Partition::PrintTables(bool print_bad) {
  // Print Tables.
  std::cout << "Good: " << std::endl;

  auto lt = good_.lock_table();
  for (auto &item : lt)
    std::cout << "Log Likelihood for " << item.first << " : " << item.second
              << std::endl;

  if (print_bad) {
    std::cout << "Bad: " << std::endl;

    auto lt1 = bad_.lock_table();
    for (auto &item : lt1)
      std::cout << "Log Likelihood for " << item.first << " : " << item.second
                << std::endl;
  }
}
/// @brief Perform both NNI moves at an edge and compare their log-likelihoods
/// to the ML, sort accordingly.
/// @param[in] tree
/// The current edge.
/// @param[in] lambda
/// The likelihood of ML tree.
/// @param[in] cutoff
/// The scaler cutoff for tree acceptance (c * lambda)
void Partition::NNIComputeEdge(pll_utree_t *tree, double lambda,
                               double cutoff) {
  // Create a clone of the original tree to perform NNI and reordering on.
  pll_utree_t *clone = pll_utree_clone(tree);

  // Set scaler parameter to determine if tree is good/bad.
  double c = cutoff;
  // Perform first NNI and reordering on first edge.
  clone = NNIUpdate(clone, 1);
  std::string label = ToNewick(clone);
  if (!(good_.contains(label) || bad_.contains(label))) {
    /// FullBranchOpt(clone);
    double lambda_1 = FullTraversalLogLikelihood(clone);
    // Compare new likelihood to ML, then decide which table to put in.
    if (lambda_1 > c * lambda) {
      good_.insert(ToNewick(clone), lambda_1);
      //Create thread for good tree and have it MakeTables.
      vec_thread_.push_back(
          std::thread(&pt::Partition::MakeTables, this, c, lambda, clone));
      //Let threads run in paralel.
      vec_thread_.at(vec_thread_.size() - 1).detach();
    } else {
      bad_.insert(ToNewick(clone), lambda_1);
    }
  }
  // Repeat for 2nd NNI move.
  clone = pll_utree_clone(tree);
  clone = NNIUpdate(clone, 2);
  label = ToNewick(clone);
  if (!(good_.contains(label) || bad_.contains(label))) {
    /// FullBranchOpt(clone);
    double lambda_1 = FullTraversalLogLikelihood(clone);
    if (lambda_1 > c * lambda) {
      good_.insert(ToNewick(clone), lambda_1);
      //Create thread for good tree and have it MakeTables.
      vec_thread_.push_back(
          std::thread(&pt::Partition::MakeTables, this, c, lambda, clone));
      //Let threads run in paralel.
      vec_thread_.at(vec_thread_.size() - 1).detach();
    } else {
      bad_.insert(ToNewick(clone), lambda_1);
    }
  }
}
/// @brief Traverse the tree and perform NNI moves at each internal edge.
void Partition::NNITraverse(pll_utree_t *tree, double lambda, double cutoff) {
  if (!tree->next)
    return;
  NNIComputeEdge(tree, lambda, cutoff);
  NNITraverse(tree->next->back, lambda, cutoff);
  NNITraverse(tree->next->next->back, lambda, cutoff);
}
}
