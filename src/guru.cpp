#include "guru.hpp"

#include <string>
#include <vector>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

#include <model_parameters.hpp>
#include <pll_partition.hpp>
#include <pll_util.hpp>

#include "wanderer.hpp"

namespace pt {

Guru::Guru(double lnl_offset,
           size_t thread_count,
           pll_utree_t* starting_tree,
           unsigned int tip_node_count,
           const pll::ModelParameters& model_parameters,
           const std::vector<std::string>& labels,
           const std::vector<std::string>& sequences,
           bool try_all_moves) :
    authority_(0.0, lnl_offset),
    thread_count_(thread_count),
    tip_node_count_(tip_node_count),
    model_parameters_(model_parameters),
    labels_(labels),
    sequences_(sequences),
    try_all_moves_(try_all_moves),
    partition_(starting_tree, tip_node_count, model_parameters, labels, sequences)
{
  AddStartingTree(starting_tree);

  // use the starting tree's log-likelihood as the authority's initial maximum
  partition_.TraversalUpdate(starting_trees_.front(), pll::TraversalType::FULL);
  authority_.SetMaximum(partition_.LogLikelihood(starting_trees_.front()));
}

Guru::~Guru()
{
  // TODO: if starting_trees_ isn't empty when the guru is destroyed,
  //       something went wrong, but it's not kosher to throw an
  //       exception from a destructor. what should we do?

  while (!starting_trees_.empty()) {
    pll_utree_every(starting_trees_.front(), pll::cb_erase_data);
    pll_utree_destroy(starting_trees_.front());
    starting_trees_.pop();
  }
}

void Guru::AddStartingTree(pll_utree_t* starting_tree)
{
  // we don't want to take ownership of starting_tree, so clone it
  // first and push the clone onto the queue
  pll_utree_t* tree = pll_utree_clone(starting_tree);
  pll_utree_every(tree, pll::cb_copy_clv_traversal);

  starting_trees_.push(tree);
}

void Guru::Start()
{
  while (wanderers_.size() < thread_count_ && !starting_trees_.empty()) {
    pll_utree_t* tree = starting_trees_.front();

    wanderers_.emplace_back(authority_, tree, tip_node_count_,
                            model_parameters_, labels_, sequences_,
                            try_all_moves_);

    // the wanderer will clone the tree for itself, so we're done with it
    pll_utree_every(tree, pll::cb_erase_data);
    pll_utree_destroy(tree);

    starting_trees_.pop();
  }

  // TODO: parallelize this
  for (auto& wanderer : wanderers_) {
    wanderer.Start();
  }
}

Authority& Guru::GetAuthority()
{
  return authority_;
}

} // namespace pt
