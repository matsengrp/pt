#include "guru.hpp"

#include <future>
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
    Authority(0.0, lnl_offset),
    thread_count_(thread_count),
    tip_node_count_(tip_node_count),
    model_parameters_(model_parameters),
    labels_(labels),
    sequences_(sequences),
    try_all_moves_(try_all_moves),
    partition_(starting_tree, tip_node_count, model_parameters, labels, sequences),
    default_tree_(nullptr)
{
  // when the guru starts, if we don't have enough starting trees to
  // create as many wanderers as requested, we'll need a default tree
  // to initialize them with
  default_tree_ = pll_utree_clone(starting_tree);
  pll_utree_every(default_tree_, pll::cb_copy_clv_traversal);

  // use the default tree's log-likelihood as the authority's initial maximum
  partition_.TraversalUpdate(default_tree_, pll::TraversalType::FULL);
  SetMaximumScore(partition_.LogLikelihood(default_tree_));

  // add the starting tree to the queue
  AddStartingTree(starting_tree);
}

Guru::~Guru()
{
  pll_utree_every(default_tree_, pll::cb_erase_data);
  pll_utree_destroy(default_tree_);

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

    wanderers_.emplace_back(*this, tree, tip_node_count_,
                            model_parameters_, labels_, sequences_,
                            try_all_moves_);

    // if an earlier wanderer happens to move to this wanderer's
    // starting tree before this wanderer is started, the wanderer's
    // Start() method will return immediately and it will go idle.

    auto& wanderer = wanderers_.back();
    futures_.emplace_back(std::async(std::launch::async,
                                     [&wanderer]() { wanderer.Start(); }));

    // the wanderer will clone the tree for itself, so we're done with it
    pll_utree_every(tree, pll::cb_erase_data);
    pll_utree_destroy(tree);

    starting_trees_.pop();
  }

  // initialize the wanderers that will start idle at the default tree
  for (size_t i = wanderers_.size(); i < thread_count_; ++i)
  {
    wanderers_.emplace_back(*this, default_tree_, tip_node_count_,
                            model_parameters_, labels_, sequences_,
                            try_all_moves_);

    // we "start" the wanderer here only to create its future. Start()
    // will immediately return when it sees that its starting tree has
    // already been visited, so its future will almost immediately be
    // in "ready" status.
    auto& wanderer = wanderers_.back();
    futures_.emplace_back(std::async(std::launch::async,
                                     [&wanderer]() { wanderer.Start(); }));
  }
}

void Guru::Wait()
{
  for (auto& future : futures_) {
    future.wait();

    // ensure exceptions are propagated
    future.get();
  }
}

} // namespace pt
