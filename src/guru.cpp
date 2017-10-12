#include "guru.hpp"

#include <algorithm>
#include <chrono>
#include <future>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <libpll/pll.h>
#include <model.hpp>
#include <pll_partition.hpp>
#include <pll_util.hpp>

#include "common.hpp"
#include "compressed_tree.hpp"
#include "move_tester.hpp"
#include "options.hpp"
#include "ordered_tree.hpp"
#include "wanderer.hpp"

// TODO: debugging only
#include <iostream>

namespace pt {

Guru::Guru(const Options& options,
           pll_utree_t* starting_tree,
           const pll::Model& model,
           const std::vector<std::string>& labels,
           const std::vector<std::string>& sequences) :
    Authority(options, 0.0),
    thread_count_(options.thread_count),
    model_(model),
    labels_(labels),
    sequences_(sequences),
    move_tester_(options.move_tester),
    partition_(starting_tree, model, labels, sequences),
    default_tree_(nullptr),
    idle_wanderer_count_(0),
    wanderer_ready_(options.thread_count, false)
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
  // synchronized in AddUnsafeStartingTree(). see comments in
  // AddUnsafeStartingTree().
  AddSafeStartingTree(default_tree_);
}

Guru::Guru(const Options& options,
           const std::vector<pll_utree_t*>& starting_trees,
           const pll::Model& model,
           const std::vector<std::string>& labels,
           const std::vector<std::string>& sequences) :
    Authority(options, 0.0),
    thread_count_(options.thread_count),
    model_(model),
    labels_(labels),
    sequences_(sequences),
    move_tester_(options.move_tester),
    partition_(starting_trees.at(0), model, labels, sequences),
    default_tree_(nullptr),
    idle_wanderer_count_(0),
    wanderer_ready_(options.thread_count, false)
{
  //
  // it is imperative that the tip node partition data indices of the
  // starting trees are synchronized to the first tree in the vector
  // before they are passed to the guru constructor. otherwise, the
  // behavior of the partition when operating on trees other than the
  // first is undefined, since the partition maintains state based on
  // the tips of the tree it was constructed with.
  //

  // clone the starting trees, because SortStartingTrees() performs a
  // traversal on each tree that will allocate node data for each node
  // if it's not already allocated, and we don't want to change the
  // trees we were given. we'll just add the clones directly to the
  // starting tree queue instead of using AddSafeStartingTree(), so
  // there's no additional cost to doing it this way.
  std::vector<pll_utree_t*> starting_clones;
  for (auto tree : starting_trees) {
    pll_utree_t* clone = pll_utree_clone(tree);
    pll_utree_every(clone, pll::cb_copy_clv_traversal);

    starting_clones.push_back(clone);
  }

  // sort the starting tree clones by descending log-likelihood
  std::vector<pll_utree_t*> sorted_clones = SortStartingTrees(starting_clones);

  // when the guru starts, if we don't have enough starting trees to
  // create as many wanderers as requested, we'll need a default tree
  // to initialize them with
  default_tree_ = pll_utree_clone(sorted_clones[0]);
  pll_utree_every(default_tree_, pll::cb_copy_clv_traversal);

  // add the starting trees to the queue
  for (auto clone : sorted_clones) {
    // since we own these clones, we can just add them to the queue
    // directly instead of using AddSafeStartingTree()
    starting_trees_.push(clone);
  }

  // use the default tree's log-likelihood as the authority's initial maximum
  pll_unode_t* root = GetVirtualRoot(default_tree_);
  partition_.TraversalUpdate(root, pll::TraversalType::FULL);
  SetMaximumScore(partition_.LogLikelihood(root));
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

  {
    std::lock_guard<std::mutex> lock(mutex_);

    starting_trees_.push(tree);
  }
}

void Guru::AddUnsafeStartingTree(pll_utree_t* starting_tree)
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

  {
    std::lock_guard<std::mutex> lock(mutex_);

    starting_trees_.push(tree);
  }
}

std::vector<pll_utree_t*> Guru::SortStartingTrees(
    const std::vector<pll_utree_t*>& starting_trees)
{
  using ValueType = std::pair<pll_utree_t*, double>;
  std::vector<ValueType> tree_lnls;

  for (auto tree : starting_trees) {
    pll_unode_t* root = GetVirtualRoot(tree);
    partition_.TraversalUpdate(root, pll::TraversalType::FULL);
    double lnl = partition_.LogLikelihood(root);

    tree_lnls.emplace_back(tree, lnl);
  }

  std::sort(tree_lnls.begin(), tree_lnls.end(),
            [](const ValueType& lhs, const ValueType& rhs) {
              return lhs.second > rhs.second;
            });

  std::vector<pll_utree_t*> sorted_trees;
  for (auto value : tree_lnls) {
    sorted_trees.emplace_back(value.first);
  }

  return sorted_trees;
}

bool Guru::RequestMove(pll_utree_t* tree, pll_unode_t* node, MoveType type)
{
  // apply the move
  pll_utree_nni(node, type, nullptr);

  bool request_accepted;

  if (idle_wanderer_count_ > 0) {
    // steal the tree. we know this tree is safe (i.e. its tip
    // node/partition data indices are synchronized with the default
    // tree and thus the wanderer partitions) because the wanderer
    // must have arrived at this tree along a path of safe trees.
    AddSafeStartingTree(tree);
    request_accepted = false;
  } else {
    CompressedTree key = GetKey(tree);
    request_accepted = GetVisitedTreeTable().insert(key, 0.0);
  }

  // undo the move
  pll_utree_nni(node, type, nullptr);

  return request_accepted;
}

void Guru::Start()
{
  // prevent launched wanderers from entering any other critical
  // section in the guru until we're done launching all of them
  std::lock_guard<std::mutex> lock(mutex_);

  while (wanderers_.size() < thread_count_ && !starting_trees_.empty()) {
    pll_utree_t* tree = starting_trees_.front();
    starting_trees_.pop();

    wanderers_.emplace_back(*this, tree,
                            model_, labels_, sequences_,
                            move_tester_);

    // if an earlier wanderer happens to move to this wanderer's
    // starting tree before this wanderer is started, the wanderer's
    // Start() method will return immediately and it will go idle.

    auto& wanderer = wanderers_.back();
    futures_.emplace_back(std::async(std::launch::async,
                                     [&wanderer]() { wanderer.Start(); }));

    // the wanderer will clone the tree for itself, so we're done with it
    pll_utree_destroy(tree, pll::cb_erase_data);
  }

  // initialize the wanderers that will start idle at the default tree
  for (size_t i = wanderers_.size(); i < thread_count_; ++i)
  {
    wanderers_.emplace_back(*this, default_tree_,
                            model_, labels_, sequences_,
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

      pll_utree_t* tree = nullptr;
      {
        std::lock_guard<std::mutex> lock(mutex_);

        if (starting_trees_.empty()) {
          continue;
        }

        tree = starting_trees_.front();
        starting_trees_.pop();
      }

      //
      // we have a tree for the idle wanderer to work on
      //

      auto& wanderer = wanderers_[i];

      wanderer.Teleport(tree);
      futures_[i] = std::async(std::launch::async,
                               [&wanderer]() { wanderer.Start(); });

      wanderer_ready_[i] = false;
      --idle_wanderer_count_;

      //std::cerr << "wanderer " << i << " launched on " << GetKey(tree) << "\n";

      // the wanderer will clone the tree for itself, so we're done with it
      pll_utree_destroy(tree, pll::cb_erase_data);
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(100));
  }

  if (!starting_trees_.empty()) {
    throw std::runtime_error("starting_trees_ is not empty");
  }
}

} // namespace pt
