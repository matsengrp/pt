#include "wanderer.hpp"

#include <queue>
#include <vector>

namespace pt {

//
// Authority
//

Authority::Authority(double ml_lnl, double lambda) :
    ml_lnl_(ml_lnl),
    lambda_(lambda)
{ }

double Authority::GetThreshold() const
{
  return ml_lnl_ * lambda_;
}

void Authority::SetThreshold(double lnl)
{
  ml_lnl_ = lnl;
}

//
// Wanderer
//

Wanderer::Wanderer(const Authority& authority, pll_partition_t* partition,
                   pll_utree_t* tree) :
    authority_(authority),
    partition_(partition),
    tree_(tree)
{ }

Wanderer::~Wanderer()
{
  pll_utree_every(tree_, cb_erase_data);
  pll_utree_destroy(tree_);

  pll_partition_destroy(partition_);
}

void Wanderer::Start()
{
  UpdatePartition(tree_, TraversalType::FULL);
  QueueMoves();

  while (move_queues_.size() > 0) {
    while (move_queues_.top().size() > 0) {
      TreeMove next_move = move_queues_.top().pop();
      Move(next_move);
    }

    // we have an empty move queue for the current tree
    move_queues_.pop();
    MoveBack();
  }

}

bool Wanderer::Move(const TreeMove& move)
{
  // do we want to store the move here, or wait until it's actually
  // accepted later?
  move_history_.push(move);
  pll_utree_nni(move.node_, move.type_, nullptr);

  // unlike the existing ToOrderedNewick() function in pll-utils.hpp,
  // this one should clone the tree, reorder it lexicographically, and
  // return a newick string (after freeing the clone). reordering
  // tree_ in-place runs the risk of making the rollback move invalid,
  // since the branches leading away from the edge may not be in the
  // same places as they were when the move was made.
  std::string newick_str = ToOrderedNewick(tree_);

  if (all_trees_.contains(newick_str)) {
    // as below it might not make sense to use MoveBack() here, since
    // we just need to reverse one move and restore a branch length,
    // then bring the partition back in sync. maybe this function
    // should be called "TryMove()" and we could have a "RejectMove()"
    // method that resets it back to the way it was if TryMove()
    // returns false?
    MoveBack();
    return false;
  }

  // invalidate CLVs around edge
  // ...

  double test_lnl = OptimizeBranch(move.node_);
  // how do we restore this branch length when we move back? should
  // OptimizeBranch() also ensure that the CLVs are oriented back
  // toward our "root"?

  all_trees_.insert(newick_str, test_lnl);

  if (test_lnl < authority_.GetThreshold()) {
    // it might make sense not to use MoveBack() here, since we really
    // only need to reverse a single move and restore a single branch
    // length, then bring the partition back in sync with a partial
    // traversal. maybe MoveBack() could be reserved for when a move
    // is actually accepted and further moves from that tree are
    // queued up.

    MoveBack();
    return false;
  }

  double lnl = OptimizeAllBranches();
  // how do we restore all the branch lengths when we move back?
  // should OptimizeAllBranches() also ensure that the CLVs are
  // oriented back toward our "root"?

  // if we find a new maximum likelihood, update the authority. note
  // there could be a race between getting and setting the maximum
  // here. maybe we need an Authority::SetMaximumIfBetter() with a
  // proper mutex or something?
  if (lnl > authority_.GetMaximum()) {
    authority_.SetMaximum(lnl);
  }

  good_trees_.insert(newick_str, lnl);
  QueueMoves();

  return true;
}

void Wanderer::MoveBack()
{
  // this function must perfectly restore the tree and partition to
  // its previous state, so it must ensure that branch lengths are all
  // restored as well as rolling back the NNI move

  // would it just make more sense to clone the tree before a move?
  // we'd have to worry about getting the partition back in sync again
  // after popping a tree off the stack, which gets us back to the
  // problem of requiring a full traversal every time a move (or step
  // back) is made. though I think technically we'd only have to
  // update the partition again if there are moves left to make on the
  // restored tree, otherwise we just pop that tree off the stack too.
}

void Wanderer::QueueMoves()
{
  std::vector<pll_utree_t*> inner_nodes(partition_->tips - 2, nullptr);
  pll_utree_query_innernodes(tree_, inner_nodes.data());

  std::queue<TreeMove> available_moves;
  for (auto node : inner_nodes) {
    // skip any pendant edges
    if (!node->back->next) {
      continue;
    }

    available_moves.push(TreeMove{node, MoveType::LEFT});
    available_moves.push(TreeMove{node, MoveType::RIGHT});
  }

  move_queues_.push(std::move(available_moves));
}

// this function should only be called when move_queues_ and
// move_history_ are empty, i.e. the wanderer is idle
void Wanderer::Teleport(pll_utree_t* tree)
{
  pll_utree_every(tree_, cb_erase_data);
  pll_utree_destroy(tree_);

  tree_ = pll_utree_clone(tree);

  Start();
}

} // namespace pt
