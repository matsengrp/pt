#include "pll_partition.hpp"

#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

#include "model_parameters.hpp"
#include "pll-utils.hpp"

namespace pt { namespace pll {

Partition::Partition(pll_utree_t* tree, unsigned int tip_node_count,
                     const ModelParameters& parameters,
                     const std::vector<std::string>& labels,
                     const std::vector<std::string>& sequences) :
    tip_node_count_(tip_node_count),
    partition_(nullptr, &pll_partition_destroy),
    sumtable_(nullptr)
{
  // TODO: accepting tip_node_count as a constructor parameter is
  //       redundant, since it's also the number of labels and
  //       sequences. might be useful as a sanity check, but the
  //       caller could be in charge of that...

  // FIXME: the variable tip_node_count shadows the member function
  //        tip_node_count() here. I'm okay with this for now as the
  //        variable is not long for this world.

  // tree is already parsed

  // check that tree is healthy (i.e., every branch has a length)
  //
  // TODO: should the caller be responsible for checking the tree?
  if (!TreeHealthy(tree)) {
    throw std::invalid_argument("Tree is missing branch lengths");
  }

  const unsigned int site_count = sequences[0].size();

  // TODO: replace macros/global variables with enums or arguments as
  //       appropriate. should they go in the model parameters?
  partition_ = PartitionPtr(
      pll_partition_create(tip_node_count_,     // tips
                           inner_node_count(),  // clv_buffers
                           STATES,              // states
                           site_count,          // sites
                           1,                   // rate_matrices
                           branch_count(),      // prob_matrices
                           RATE_CATS,           // rate_cats
                           inner_node_count(),  // scale_buffers
                           ARCH_FLAGS),         // attributes
      &pll_partition_destroy);

  // TODO: try/catch?
  // TODO: should these be static member functions? does it matter?
  SetModelParameters(parameters);
  SetTipStates(tree, labels, sequences);

  // TODO: is there a better place for this? as far as I can tell from
  //       the docs, this array is never updated, so it could probably
  //       be a const vector and initialized in the initializer list
  params_indices_.assign(RATE_CATS, 0);

  AllocateScratchBuffers();
}

Partition::Partition(Partition&& rhs) :
    tip_node_count_(rhs.tip_node_count_),
    partition_(nullptr, &pll_partition_destroy),
    sumtable_(nullptr)
{
  partition_ = std::move(rhs.partition_);
  params_indices_.assign(RATE_CATS, 0);

  AllocateScratchBuffers();
}

Partition::~Partition()
{
  FreeScratchBuffers();
}

void Partition::SetModelParameters(const ModelParameters& parameters)
{
  std::vector<double> rate_cats(RATE_CATS, 0.0);

  pll_compute_gamma_cats(parameters.alpha,
                         rate_cats.size(),
                         rate_cats.data());

  // set frequencies at model with index 0 (we currently have only one model).
  pll_set_frequencies(partition_.get(), 0, parameters.frequencies.data());

  // set 6 substitution parameters at model with index 0
  pll_set_subst_params(partition_.get(), 0, parameters.subst_params.data());

  // set rate categories
  pll_set_category_rates(partition_.get(), rate_cats.data());
}

void Partition::SetTipStates(pll_utree_t* tree,
                             const std::vector<std::string>& labels,
                             const std::vector<std::string>& sequences)
{
  if (labels.size() != sequences.size()) {
    throw std::invalid_argument("Number of labels does not match number of sequences");
  }

  // TODO: if tip_node_count_ is inferred from the size of the labels
  //       or sequences vector in the constructor, this check becomes
  //       redundant
  if (tip_node_count_ != labels.size()) {
    throw std::invalid_argument("Unexpected number of tip nodes supplied");
  }

  // obtain an array of pointers to tip nodes
  std::vector<pll_utree_t*> tip_nodes(tip_node_count_, nullptr);
  pll_utree_query_tipnodes(tree, tip_nodes.data());

  std::map<std::string, unsigned int> tip_ids;

  // populate a hash table with tree tip labels
  for (unsigned int i = 0; i < tip_node_count_; ++i) {
    std::string label = tip_nodes[i]->label;
    bool inserted;
    std::tie(std::ignore, inserted) = tip_ids.emplace(label, i);

    if (!inserted) {
      throw std::invalid_argument("Error inserting tip label " + label
                                  + " into map (possibly a duplicate)");
    }
  }

  // find sequences in hash table and link them with the corresponding taxa
  for (unsigned int i = 0; i < tip_node_count_; ++i) {
    auto iter = tip_ids.find(labels[i]);

    if (iter == tip_ids.end()) {
      throw std::invalid_argument("Sequence with header " + labels[i]
                                  + " does not appear in the tree");
    }

    unsigned int tip_clv_index = iter->second;
    pll_set_tip_states(partition_.get(), tip_clv_index, pll_map_nt, sequences[i].c_str());
  }
}

void Partition::AllocateScratchBuffers()
{
  // allocate scratch buffers for TraversalUpdate()
  travbuffer_.resize(node_count());
  branch_lengths_.resize(branch_count());
  matrix_indices_.resize(branch_count());
  operations_.resize(inner_node_count());

  // allocate scratch buffers for OptimizeBranch()
  sumtable_ = (double*) pll_aligned_alloc(
      partition_->sites * partition_->rate_cats * partition_->states_padded * sizeof(double),
      ALIGNMENT);
}

void Partition::FreeScratchBuffers()
{
  if (sumtable_) {
    pll_aligned_free(sumtable_);
    sumtable_ = nullptr;
  }
}

double Partition::LogLikelihood(pll_utree_t* tree)
{
  double lnl = pll_compute_edge_loglikelihood(
      partition_.get(), tree->clv_index, tree->scaler_index,
      tree->back->clv_index, tree->back->scaler_index, tree->pmatrix_index,
      params_indices_.data(), nullptr);

  return lnl;
}

unsigned int Partition::TraversalUpdate(pll_utree_t* root, TraversalType type)
{
  unsigned int traversal_size;
  int status;

  if (type == TraversalType::FULL) {
    status = pll_utree_traverse(root, cb_full_traversal, travbuffer_.data(),
                                &traversal_size);
  } else if (type == TraversalType::PARTIAL) {
    status = pll_utree_traverse(root, cb_partial_traversal, travbuffer_.data(),
                                &traversal_size);
  } else {
    throw std::invalid_argument("Invalid traversal type");
  }

  if (!status) {
    throw std::invalid_argument("TraversalUpdate() requires an inner node");
  }

  unsigned int matrix_count;
  unsigned int ops_count;

  // Given the computed traversal descriptor, generate the operations
  // structure, and the corresponding probability matrix indices that
  // may need recomputing.
  pll_utree_create_operations(travbuffer_.data(), traversal_size,
                              branch_lengths_.data(), matrix_indices_.data(),
                              operations_.data(), &matrix_count, &ops_count);

  // Update matrix_count probability matrices for model with index 0. The i-th
  // matrix (i ranges from 0 to matrix_count - 1) is generated using branch
  // length branch_lengths[i] and can be referred to with index
  // matrix_indices[i].
  pll_update_prob_matrices(partition_.get(), params_indices_.data(),
                           matrix_indices_.data(), branch_lengths_.data(),
                           matrix_count);

  // Use the operations array to compute all ops_count inner CLVs. Operations
  // will be carried out sequentially starting from operation 0 towards
  // ops_count-1.
  pll_update_partials(partition_.get(), operations_.data(), ops_count);

  return ops_count;
}

void Partition::UpdateBranchLength(pll_utree_t* node, double length)
{
  // Update current branch lengths.
  node->length = length;
  node->back->length = length;

  // Update this branch's probability matrix now that the branch
  // length has changed. No CLVs need to be invalidated.
  pll_update_prob_matrices(partition_.get(), params_indices_.data(),
                           &(node->pmatrix_index),
                           &(node->length), 1);
}

double Partition::OptimizeBranch(pll_utree_t* node)
{
  pll_utree_t *parent = node;
  pll_utree_t *child = node->back;

  // Compute the sumtable for the particular branch once before proceeding with
  // the optimization.
  pll_update_sumtable(partition_.get(), parent->clv_index, child->clv_index,
                      params_indices_.data(), sumtable_);

  double len = node->length;
  bool maybe_decreasing = false;

  for (unsigned int i = 0; i < MAX_ITER; ++i) {
    double d1; // First derivative.
    double d2; // Second derivative.

    pll_compute_likelihood_derivatives(
        partition_.get(), parent->scaler_index, child->scaler_index, len,
        params_indices_.data(), sumtable_, &d1, &d2);

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
    if (d2 < 0.0)
      len += d1 / d2;
    else
      len -= d1 / d2;

    // If the next branch length to evaluate goes negative, we instead
    // set it to a small positive value for the next iteration. If
    // this has happened before, we stop early, as the curve is
    // probably decreasing.
    if (len < 0.0) {
      len = EPSILON;

      if (maybe_decreasing) {
        break;
      }

      maybe_decreasing = true;
    }
  }

  // No CLVs need to be invalidated; see the definition of UpdateBranchLength().
  UpdateBranchLength(node, len);

  return len;
}

void Partition::OptimizeAllBranchesOnce(pll_utree_t* tree)
{
  std::vector<pll_utree_t*> nodes(node_count(), nullptr);
  unsigned int nodes_found;

  // Traverse the entire tree and collect nodes using a callback
  // function that returns 1 for every node visited. Some of these
  // nodes will be tips, in which case we operate on node->back (the
  // tip's parent) instead of node; see below.
  if (!pll_utree_traverse(tree, [](pll_utree_t*) { return 1; }, nodes.data(),
                          &nodes_found)) {
    throw std::invalid_argument("OptimizeAllBranches() requires an inner node");
  }

  if (nodes_found != nodes.size()) {
    throw std::invalid_argument("Unexpected number of nodes");
  }

  for (auto node : nodes) {
    // If this is a tip node, operate on its parent instead.
    if (!node->next) {
      node = node->back;
    }

    TraversalUpdate(node, TraversalType::PARTIAL);
    OptimizeBranch(node);
  }

  // Leave the tree and partition in the same state we found it.
  //
  // TODO: This assumes that the CLVs were originally oriented toward
  // the node pointed to by tree -- this may not be the case, or the
  // caller may not care. Maybe we should just place the
  // responsibility on the caller?
  TraversalUpdate(tree, TraversalType::PARTIAL);
}

void Partition::OptimizeAllBranches(pll_utree_t* tree)
{
  // TODO: Why are we doing a full traversal here instead of partial?
  TraversalUpdate(tree, TraversalType::FULL);
  double loglike_prev = LogLikelihood(tree);

  OptimizeAllBranchesOnce(tree);

  // TODO: If we remove the call to TraversalUpdate() at the end of
  //       OptimizeAllBranchesOnce(), we'll need to do a partial
  //       traversal before computing the log-likelihood.
  double loglike = LogLikelihood(tree);

  unsigned int i = 0;
  while (fabs(loglike_prev - loglike) > EPSILON && i < MAX_ITER) {
    OptimizeAllBranchesOnce(tree);

    loglike_prev = loglike;

    // TODO: As above, if we remove the call to TraversalUpdate() at
    //       the end of OptimizeAllBranchesOnce(), we'll need to do a
    //       partial traversal before computing the log-likelihood.
    loglike = LogLikelihood(tree);
    i++;
  }
}

} } // namespace pt::pll
