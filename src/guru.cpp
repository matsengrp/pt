#include "guru.hpp"

#include <chrono>
#include <future>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <libpll/pll.h>
#include <model_parameters.hpp>
#include <pll_partition.hpp>
#include <pll_util.hpp>

#include "common.hpp"
#include "move_tester.hpp"
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
           std::shared_ptr<const MoveTester> move_tester) :
    Authority(0.0, lnl_offset),
    thread_count_(thread_count),
    model_parameters_(model_parameters),
    labels_(labels),
    sequences_(sequences),
    move_tester_(move_tester),
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
  //
  // TODO: the partition is only ever used here -- it doesn't actually
  //       need to be a member variable if that's the case
  pll_unode_t* root = GetVirtualRoot(default_tree_);
  partition_.TraversalUpdate(root, pll::TraversalType::FULL);
  SetMaximumScore(partition_.LogLikelihood(root));

  // add the starting tree to the queue. default_tree_ will be the
  // tree to which all tip node/partition data indices are
  // synchronized in the public AddStartingTree(). see comments in
  // AddStartingTree().
  AddSafeStartingTree(default_tree_);
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

void Guru::AddSafeStartingTree(pll_utree_t* starting_tree)
{
  // we don't want to take ownership of starting_tree, so clone it
  // first and push the clone onto the queue
  pll_utree_t* tree = pll_utree_clone(starting_tree);
  pll_utree_every(tree, pll::cb_copy_clv_traversal);

  starting_trees_.push(tree);
}

void Guru::AddStartingTree(pll_utree_t* starting_tree)
{
  //
  // we don't know where this tree came from; for example, if it was
  // parsed by libpll from a Newick string, there's no guarantee that
  // its tips have the same node/partition data indices as the default
  // tree the guru was constructed with. wanderer partitions maintain
  // state related to the tips, which means we need to synchronize
  // this tree's tips to those of the default tree based on the tip
  // labels before it can be considered safe for a wanderer to use.
  //

  // we don't want to take ownership of starting_tree, so clone it
  // first, synchronize the tips with the default tree, and push the
  // clone onto the queue
  pll_utree_t* tree = pll_utree_clone(starting_tree);
  pll_utree_every(tree, pll::cb_copy_clv_traversal);

  pll::SynchronizeTipIndices(default_tree_, tree);

  starting_trees_.push(tree);
}

bool Guru::RequestMove(pll_utree_t* tree, pll_unode_t* node, MoveType type)
{
  // apply the move
  pll_utree_nni(node, type, nullptr);

  std::string newick_str = GetKey(tree);
  bool request_accepted;

  if (idle_wanderer_count_ > 0) {
    std::lock_guard<std::mutex> lock(mutex_);

    // steal the tree. we know this tree is safe (i.e. its tip
    // node/partition data indices are synchronized with the default
    // tree and thus the wanderer partitions) because the wanderer
    // must have arrived at this tree along a path of safe trees.
    AddSafeStartingTree(tree);
    request_accepted = false;
  } else {
    request_accepted = GetVisitedTreeTable().insert(newick_str, 0.0);
  }

  // undo the move
  pll_utree_nni(node, type, nullptr);

  return request_accepted;
}

void Guru::Start()
{
  while (wanderers_.size() < thread_count_ && !starting_trees_.empty()) {
    pll_utree_t* tree = starting_trees_.front();

    wanderers_.emplace_back(*this, tree,
                            model_parameters_, labels_, sequences_,
                            move_tester_);

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
                            move_tester_);

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
