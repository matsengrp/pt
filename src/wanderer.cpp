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
    current_tree_(tree)
{ }

Wanderer::~Wanderer()
{
  pll_utree_every(current_tree_, cb_erase_data);
  pll_utree_destroy(current_tree_);

  pll_partition_destroy(partition_);
}

void Wanderer::Start()
{
  UpdatePartition(current_tree_, TraversalType::FULL);
  QueueMoves();

  // note that move_queues_.top() and move_queues_.top().front() can
  // change during iteration, so we don't store references to them
  // inside the loops. we always want to be looking at the next move
  // in the top queue.
  while (move_queues_.size() > 0) {
    while (move_queues_.top().size() > 0) {
      // somewhere here is where we would check to see if there are
      // any idle monks, and use Teleport() (or ask the authority to)
      // to have the idle monk go there instead of taking the move
      // ourselves

      MoveForward();
      QueueMoves();
    }

    // we have an empty move queue for the current tree
    move_queues_.pop();
    MoveBack();
  }

}

bool Wanderer::TryMove(pll_utree_t* node, MoveType type)
{
  pll_utree_nni(node, type, nullptr);

  std::string newick_str = ToOrderedNewick(current_tree_);

  // it would probably be better here if this was an insert rather
  // than a contains, so that we know for certain that this Wanderer
  // "owns" the tree it's currently evaluating. otherwise the code
  // that follows may be evaluated by multiple monks
  if (all_trees.contains(newick_str)) {
    // undo the move and reject
    pll_utree_nni(node, type, nullptr);
    return false;
  }

  // if PT_TEST_SINGLE_BRANCH is defined, TryMove() will optimize the
  // branch length and use the resulting log-likelihood to determine
  // if this move should be accepted. if it's not defined, TryMove()
  // will always return true if this tree hasn't been visited before.
#ifdef PT_TEST_SINGLE_BRANCH
  //
  // so the big question here is, when do we use current_tree_ as the argument
  // to the various functions below, and when do we use node?
  //

  const double original_length = node->length;

  OptimizeBranch(node);

  // when OptimizeBranch() returns, will the CLVs be oriented toward
  // the optimized edge, or will they be pointing back to current_tree_ (our
  // "root")? what CLVs need to be invalidated, in either case? my
  // belief is that if we use the optimized edge for computing the
  // log-likelihood, no CLVs need to be invalidated. if we use the
  // "root", the CLVs around that edge must be invalidated first.

  // invalidate CLVs?
  // ...

  const double test_lnl = LogLikelihood(node);

  bool accept_move = false;
  if (test_lnl > authority.GetThreshold()) {
    accept_move = true;
  }

  // store the test log-likelihood along with the tree in the "all"
  // table; it might be useful later for comparison
  all_trees_.insert(newick_str, test_lnl);

  // restore the branch length, undo the move, invalidate any CLVs
  // (?), and do a partial traversal to restore the partition to its
  // pre-TryMove() state (?)

  UpdateBranchLength(node, original_length);
  pll_utree_nni(node, type, nullptr);

  // invalidate CLVs? maybe UpdateBranchLength() does this
  // ...

  UpdatePartition(current_tree_, TraversalType::PARTIAL); // ?

  return accept_move;
#else
  // undo the move
  pll_utree_nni(node, type, nullptr);

  // store the tree in the "all" table; we don't have any log-likelihood for it
  all_trees_.insert(newick_str, 0.0);

  return true;
#endif
}

void Wanderer::MoveForward()
{
  // clone tree onto the tree history stack

  // apply move. the move we're applying should be
  // move_queues_.top().front(), I think

  // do full branch optimization

  // compute log-likelihood

  // add tree to "good" table

  // should QueueMoves() be called in here, or in the caller?
}

void Wanderer::MoveBack()
{
  // restore previous tree by popping the current tree off the top of
  // the tree history stack

  // resynchronize partition and tree with a full traversal? this may
  // be pointless if the move queue for the previous tree is empty
}




bool Wanderer::Move(const TreeMove& move)
{
  pll_utree_nni(move.node_, move.type_, nullptr);

  // unlike the existing ToOrderedNewick() function in pll-utils.hpp,
  // this one should clone the tree, reorder it lexicographically, and
  // return a newick string (after freeing the clone). reordering
  // current_tree_ in-place runs the risk of making the rollback move invalid,
  // since the branches leading away from the edge may not be in the
  // same places as they were when the move was made.
  std::string newick_str = ToOrderedNewick(current_tree_);

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
  // is n_tips - 2 correct for the number of inner nodes? probably
  // should add an error check to see if the number of nodes
  // pll_utree_query_innernodes() finds is the same as the size of the
  // vector
  std::vector<pll_utree_t*> inner_nodes(partition_->tips - 2, nullptr);
  pll_utree_query_innernodes(current_tree_, inner_nodes.data());

  std::queue<TreeMove> available_moves;
  for (auto node : inner_nodes) {
    // skip any pendant edges
    if (!node->back->next) {
      continue;
    }

    // here is where we can try applying a move, optimizing the
    // current branch length, and then comparing the log-likelihood to
    // the threshold to see if the move is worth taking (since if we
    // do take the move, we'll be performing a full branch length
    // optimization). if the move is worth taking, we add the move to
    // the move queue. regardless, after testing the move, we undo the
    // move and restore the original branch length, then repeat for
    // the move in the other "direction".

    // assume we have a TryMove() method that applies the move,
    // optimizes the branch length, and compares the log-likelihood to
    // the authority's threshold. it then reverses the move and
    // restores the branch length. it may be okay if the CLVs are
    // oriented differently than they were before TryMove() gets
    // called, since TryMove() will have to do a partial traversal
    // anyway, right? TryMove() will return true if the move is worth
    // taking, false otherwise.

    // TryMove() could probably also be the place where we check to
    // see if the tree has already been tested/visited.

    // if the single-branch optimization behavior isn't desired,
    // TryMove() could instead just always return true, and the move
    // will be taken and full branch optimization will be performed.

    if (TryMove(node, MoveType::LEFT)) {
      available_moves.push(TreeMove{node, MoveType::LEFT});
    }

    if (TryMove(node, MoveType::RIGHT)) {
      available_moves.push(TreeMove{node, MoveType::RIGHT});
    }
  }

  // note: if we're using the TryMove() strategy described above,
  // we'll need to handle the case where the available move queue is
  // empty, both here and in the calling function.
  move_queues_.push(std::move(available_moves));
}

// this function should only be called when move_queues_ and
// move_history_ are empty, i.e. the wanderer is idle
void Wanderer::Teleport(pll_utree_t* tree)
{
  pll_utree_every(current_tree_, cb_erase_data);
  pll_utree_destroy(current_tree_);

  current_tree_ = pll_utree_clone(tree);

  Start();
}

} // namespace pt
