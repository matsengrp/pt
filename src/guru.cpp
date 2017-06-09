#include "guru.hpp"

#include <chrono>
#include <future>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <libpll/pll.h>
#include <model_parameters.hpp>
#include <pll_partition.hpp>
#include <pll_util.hpp>

#include "ordered_tree.hpp"
#include "wanderer.hpp"

// TODO: debugging only
#include <iostream>

namespace pt {

Guru::Guru(double lnl_offset,
           size_t thread_count,
           pll_utree_t* starting_tree,
           const pll::ModelParameters& model_parameters,
           const std::vector<std::string>& labels,
           const std::vector<std::string>& sequences,
           bool try_all_moves) :
    Authority(0.0, lnl_offset),
    thread_count_(thread_count),
    model_parameters_(model_parameters),
    labels_(labels),
    sequences_(sequences),
    try_all_moves_(try_all_moves),
    partition_(starting_tree, model_parameters, labels, sequences),
    default_tree_(nullptr),
    idle_wanderer_count_(0),
    wanderer_ready_(thread_count, false)
{
  // when the guru starts, if we don't have enough starting trees to
  // create as many wanderers as requested, we'll need a default tree
  // to initialize them with
  default_tree_ = pll_utree_clone(starting_tree);
  pll_utree_every(default_tree_, pll::cb_copy_clv_traversal);

  // use the default tree's log-likelihood as the authority's initial maximum
  pll_unode_t* root = GetVirtualRoot(default_tree_);
  partition_.TraversalUpdate(root, pll::TraversalType::FULL);
  SetMaximumScore(partition_.LogLikelihood(root));

  // add the starting tree to the queue
  AddStartingTree(starting_tree);
}

Guru::~Guru()
{
  pll_utree_destroy(default_tree_, pll::cb_erase_data);

  // TODO: if starting_trees_ isn't empty when the guru is destroyed,
  //       something went wrong, but it's not kosher to throw an
  //       exception from a destructor. what should we do?

  while (!starting_trees_.empty()) {
    pll_utree_destroy(starting_trees_.front(), pll::cb_erase_data);
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

std::pair<bool, std::string> Guru::RequestTree(pll_utree_t* tree, bool first_tree)
{
  std::string newick_str = GetKey(tree);

  if (GetVisitedTreeTable().contains(newick_str)) {
    return std::make_pair(false, newick_str);
  }

  // if this is the wanderer's first tree, don't steal it
  if (!first_tree && idle_wanderer_count_ > 0) {
    std::lock_guard<std::mutex> lock(mutex_);
    AddStartingTree(tree);

    return std::make_pair(false, newick_str);
  }

  return std::make_pair(GetVisitedTreeTable().insert(newick_str, 0.0), newick_str);
}

void Guru::Start()
{
  while (wanderers_.size() < thread_count_ && !starting_trees_.empty()) {
    pll_utree_t* tree = starting_trees_.front();

    wanderers_.emplace_back(*this, tree,
                            model_parameters_, labels_, sequences_,
                            try_all_moves_);

    // if an earlier wanderer happens to move to this wanderer's
    // starting tree before this wanderer is started, the wanderer's
    // Start() method will return immediately and it will go idle.

    auto& wanderer = wanderers_.back();
    futures_.emplace_back(std::async(std::launch::async,
                                     [&wanderer]() { wanderer.Start(); }));

    // the wanderer will clone the tree for itself, so we're done with it
    pll_utree_destroy(tree, pll::cb_erase_data);

    starting_trees_.pop();
  }

  // initialize the wanderers that will start idle at the default tree
  for (size_t i = wanderers_.size(); i < thread_count_; ++i)
  {
    wanderers_.emplace_back(*this, default_tree_,
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

void Guru::UpdateWandererStatus(size_t index)
{
  auto& future = futures_[index];

  if (!future.valid()) {
    // sanity check -- an invalid future indicates that we have
    // already called get() on the future, which means that we should
    // have already marked it as ready in a previous pass.
    if (!wanderer_ready_[index]) {
      throw std::runtime_error("future is invalid but wanderer not marked ready");
    }

    //std::cerr << "wanderer " << index << " still ready\n";

    return;
  }

  std::future_status status = future.wait_for(std::chrono::milliseconds(0));

  // sanity check -- no wanderer should ever be deferred
  if (status == std::future_status::deferred) {
    throw std::runtime_error("a wanderer is deferred and shouldn't be");
  }

  if (status == std::future_status::timeout) {
    wanderer_ready_[index] = false;

    //std::cerr << "wanderer " << index << " still running\n";

    return;
  }

  // sanity check
  if (status != std::future_status::ready) {
    throw std::runtime_error("unexpected future status");
  }

  //std::cerr << "wanderer " << index << " now ready\n";

  // if we're here, the wanderer went idle since it was last checked,
  // so ensure any exceptions it throws get propagated via future.get()
  // and mark it as ready for work. note that we can only ever call
  // get() on a future once, after which valid() will return false
  future.get();

  wanderer_ready_[index] = true;
  ++idle_wanderer_count_;
}

void Guru::Wait()
{
  while (idle_wanderer_count_ < thread_count_) {
    for (size_t i = 0; i < thread_count_; ++i)
    {
      UpdateWandererStatus(i);

      if (!wanderer_ready_[i]) {
        continue;
      }

      // TODO: this lock could probably be scoped more efficiently by
      //       moving the other calls to starting_trees_ methods up
      std::lock_guard<std::mutex> lock(mutex_);

      if (starting_trees_.empty()) {
        continue;
      }

      //
      // we have a tree for the idle wanderer to work on
      //

      pll_utree_t* tree = starting_trees_.front();
      auto& wanderer = wanderers_[i];

      wanderer.Teleport(tree);
      futures_[i] = std::async(std::launch::async,
                               [&wanderer]() { wanderer.Start(); });

      wanderer_ready_[i] = false;
      --idle_wanderer_count_;

      //std::cerr << "wanderer " << i << " launched on " << GetKey(tree) << "\n";

      // the wanderer will clone the tree for itself, so we're done with it
      pll_utree_destroy(tree, pll::cb_erase_data);

      starting_trees_.pop();
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }

  if (!starting_trees_.empty()) {
    throw std::runtime_error("starting_trees_ is not empty");
  }
}

} // namespace pt
