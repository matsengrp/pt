#include "partition.hpp"
#include "pll-utils.hpp"

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
Partition::Partition(std::string newick_path, std::string fasta_path, std::string RAxML_info_path) {
  // Parse the unrooted binary tree in newick format, and store the number of
  // tip nodes in tip_nodes_count.
  tree_ = pll_utree_parse_newick(&newick_path[0], &tip_nodes_count_);
  if (!tree_) {
    std::cout << "Parsing failure: " << newick_path
              << " must be an unrooted binary tree.\n";
    exit(EXIT_FAILURE);
  };

  if (!TreeHealthy(tree_)) fatal("Missing branch lengths in tree.\n");

  char **headers = NULL;
  char **seqdata = NULL;
  unsigned int sites_count =
      ParseFasta(fasta_path, tip_nodes_count(), &headers, &seqdata);

  partition_ = pll_partition_create(
      tip_nodes_count(),    // Number of tip sequences we want to have.
      inner_nodes_count(),  // Number of CLV buffers to be allocated for inner
                            // nodes.
      STATES,               // Number of states that our data have.
      sites_count,          // Number of sites in the alignment.
      1,  // Number of different substitution models (or eigen decomposition)
          //     to use concurrently (i.e. 4 for LG4).
      branch_count(),       // Number of probability matrices to be allocated.
      RATE_CATS,            // Number of rate categories we will use.
      inner_nodes_count(),  // How many scale buffers to use.
      ARCH_FLAGS);          // List of flags for hardware acceleration.

  SetModelParameters(partition_,RAxML_info_path);

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
std::string Partition::ToNewick() {
  return std::string(pll_utree_export_newick(tree_));
}


/// @brief Finds node in tree's subtree with the label that is first alphabetically.
/// @param[in] tree
/// Parent node of subtree to root.
/// @return smallest label
std::string Partition::RootNewickRecursive(pll_utree_t* tree) {
  /*Root the tree at parent of lowest-labeled leaf*/
  std::string minlabel;
  if(!tree->next)
  return tree->label;
  std::string rlabel1= RootNewickRecursive(tree->next->back);
  std::string rlabel2= RootNewickRecursive(tree->next->next->back);
  if(rlabel1.compare(rlabel2)<0)
    minlabel=rlabel1;
  else
    minlabel=rlabel2;
  return minlabel;
}


/// @brief Finds a nod with the given label and sets the current node to the parent of it.
/// @param[in] label
/// Label of node to find
/// @param[in] tree
/// Node to search at
void Partition:: SetLabelRoot(std::string label,pll_utree_t* tree){
  if(!tree->next){
    if(tree->label==label)
      tree_=tree->back;
    else
      return;
  }
  else{
    SetLabelRoot(label, tree->next->back);
    SetLabelRoot(label, tree->next->next->back);
  }
}


/// @brief Sets the label of the root for the entire tree. Sets tree_ to that node.
/// @param[in] tree
/// An internal node.
void Partition:: SetNewickRoot(pll_utree_t* tree){
  std::string minlabel;
  std::string rlabel1= RootNewickRecursive(tree);
  std::string rlabel2= RootNewickRecursive(tree->back);
  if(rlabel1.compare(rlabel2)<0)
  minlabel=rlabel1;
  else
  minlabel=rlabel2;
  SetLabelRoot(minlabel,tree);
  SetLabelRoot(minlabel,tree->back);
}


/// @brief Recursively orders all subtrees of node tree.
/// @param[in] tree
/// Node at which to recursively order.
void Partition::RecursiveOrderedNewick(pll_utree_t* tree){
  std::string newick;
  if(!tree->next)
  {
  newick.append(tree->label);
  }
  else{
    std::string tree1label;
    std::string tree2label;
    tree1label=RootNewickRecursive(tree->next->back);
    tree2label=RootNewickRecursive(tree->next->next->back);
    if(tree1label.compare(tree2label)<0){
      return;
    }
    else{
      pll_utree_t* temp= tree->next;
      tree->next=tree->next->next;
      tree->next->next=temp;
      tree->next->next->next=tree;
    }

  }
}


/// @brief Determines order for first ternary step, then runs recursive ordering for the two non-root subtrees.
/// @return Completely ordered newick string.
void Partition::ToOrderedNewick(){
  SetNewickRoot(tree_);
  std::string tree1label=RootNewickRecursive(tree_->next->back);
  std::string tree2label=RootNewickRecursive(tree_->next->next->back);
  std::string subtree1;
  std::string subtree2;
  if(tree1label.compare(tree2label)<0){
      return;
    }
  else{
  pll_utree_t* temp= tree_->next;
      tree_->next=tree_->next->next;
      tree_->next->next=temp;
      tree_->next->next->next=tree_;
  }
  /*Check tree health after reordering*/
  if(!pll_utree_check_integrity(tree_))
  fatal("Tree not healthy, check reordering parameterization.");
    if(!TreeHealthy(tree_))
  fatal("Tree not healthy, check reordering parameterization.");
}




/// @brief Perform a full tree traversal and update CLV's, etc.
/// Eventually, we will want to break this up for efficiency gains
/// so that we aren't always doing a full traversal every time we
/// want to calculate the likelihood, but we'll do that carefully!
/// @param[in] tree
/// Parent node from which to update CLV's etc.
void Partition::FullTraversalUpdate(pll_utree_t* tree) {
  /* Perform a full postorder traversal of the unrooted tree. */
  unsigned int traversal_size;
  unsigned int matrix_count, ops_count;
  if (!pll_utree_traverse(tree, cb_full_traversal, travbuffer_,
                          &traversal_size)) {
    fatal("Function pll_utree_traverse() requires inner nodes as parameters");
  }

  /* Given the computed traversal descriptor, generate the operations
     structure, and the corresponding probability matrix indices that
     may need recomputing. */
  pll_utree_create_operations(travbuffer_, traversal_size, branch_lengths_,
                              matrix_indices_, operations_, &matrix_count,
                              &ops_count);

  /* Update matrix_count probability matrices for model with index 0. The i-th
     matrix (i ranges from 0 to matrix_count - 1) is generated using branch
     length branch_lengths[i] and can be referred to with index
     matrix_indices[i]. */
  pll_update_prob_matrices(partition_, params_indices_, matrix_indices_,
                           branch_lengths_, matrix_count);

  /* Use the operations array to compute all ops_count inner CLVs. Operations
     will be carried out sequentially starting from operation 0 towrds
     ops_count-1. */
  pll_update_partials(partition_, operations_, ops_count);
}


/// @brief Just calculate log likelihood.
/// NOTE: This will return nonsense if you haven't updated CLV's, etc.
/// @return Log likelihood.
double Partition::LogLikelihood() {
  double logl = pll_compute_edge_loglikelihood(
      partition_, tree_->clv_index, tree_->scaler_index, tree_->back->clv_index,
      tree_->back->scaler_index, tree_->pmatrix_index, params_indices_,
      NULL);  // Can supply a persite_lnl parameter as last argument.

  return logl;
}


/// @brief Make a traversal at the root and return the log likelihood.
/// @return Log likelihood.
double Partition::FullTraversalLogLikelihood() {
  FullTraversalUpdate(tree_);
  return LogLikelihood();
}


/// @brief Optimize the current branch length.
/// NOTE: This assumes we've just run `FullTraversalUpdate` or some such.
/// @param[in] tree
/// Parent node of branch to optimize.
/// @return Optimized brance length.
double Partition::OptimizeCurrentBranch(pll_utree_t* tree) {
  /* Perform a full postorder traversal of the unrooted tree. */
  pll_utree_t *parent = tree;
  pll_utree_t *child = tree->back;

  // Compute the sumtable for the particular branch once before proceeding with
  // the optimization.
  pll_update_sumtable(partition_, parent->clv_index, child->clv_index,
                      params_indices_, sumtable_);

  double len = tree->length;
  double d1;  // First derivative.
  double d2;  // Second derivative.
  for (unsigned int i = 0; i < MAX_ITER; ++i) {
    double opt_logl = pll_compute_likelihood_derivatives(
        partition_, parent->scaler_index, child->scaler_index, len,
        params_indices_, sumtable_, &d1, &d2);

    printf("Branch length: %f log-L: %f Derivative: %f\n", len, opt_logl, d1);

    // If derivative is approximately zero then we've found the maximum.
    if (fabs(d1) < EPSILON) break;

    /* Newton's method for finding the optimum of a function. The iteration to
       reach the optimum is

       x_{i+1} = x_i - f'(x_i) / f''(x_i)

       where x_i is the current branch, f'(x_i) is the first derivative and
       f''(x_i) is the second derivative of the likelihood function. */
    len -= d1 / d2;

    /*Branch optimization was returning negative values, set minimum branch length to be small positive value*/

    if(len<0)
    {
    len=EPSILON;
    break;
    }
  }

  // Update current branch lengths.
  parent->length = len;
  child->length = len;

  return len;
}

/// @brief Aux function for optimizing branch lengths accross the whole tree.
/// @param[in] tree
/// Child node from which to optimize the branch length and continue recursion.
void Partition::TreeBranchLengthsAux(pll_utree_t *tree) {
  if (!tree->next){
  FullTraversalUpdate(tree->back);
  OptimizeCurrentBranch(tree->back);
  std::cout<<"DONE"<<std::endl;
  }
  else{
  FullTraversalUpdate(tree);
  OptimizeCurrentBranch(tree);
  std::cout<<"DONE"<<std::endl;

  TreeBranchLengthsAux(tree->next->back);
  TreeBranchLengthsAux(tree->next->next->back);
  }
}

///@brief Perform a postorder tree traversal, optimizing the branch length at every edge.
///NOTE: This optimizes internal branches twice per traversal.
void Partition::TreeBranchLengths(){
  if(!tree_->next)
    {
    fatal("Function TreeBranchLengthsAux requires an inner node as parameter");
    }
  TreeBranchLengthsAux(tree_);
  TreeBranchLengthsAux(tree_->back);

}

///@brief Repeat branch length optimization until log likelihood changes very little or MAX_ITER is reached.
void Partition::FullBranchOpt(){
  double loglike_prev=FullTraversalLogLikelihood();

  TreeBranchLengths();

  double loglike=FullTraversalLogLikelihood();

  int i=0;
  while((fabs(loglike_prev-loglike) > EPSILON)&& i<MAX_ITER){
  TreeBranchLengths();

  loglike_prev=loglike;
  loglike=FullTraversalLogLikelihood();
  i++;

  }

}

}
