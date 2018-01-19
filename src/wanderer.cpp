#include "wanderer.hpp"

#include <memory>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <libpll/pll.h>
#include <model.hpp>
#include <pll_partition.hpp>
#include <pll_util.hpp>

#include "authority.hpp"
#include "common.hpp"
#include "move_tester.hpp"
#include "ordered_tree.hpp"
#include "position.hpp"

namespace pt {

//
// Wanderer
//

Wanderer::Wanderer(Authority& authority,
                   const Position& starting_position,
                   const std::vector<std::string>& labels,
                   const std::vector<std::string>& sequences,
                   std::shared_ptr<const MoveTester> move_tester,
                   bool optimize_models,
                   bool map_mode) :
    authority_(authority),
    partition_(starting_position.GetTree(),
               starting_position.GetModel(),
               labels, sequences, map_mode),
    move_tester_(move_tester),
    optimize_models_(optimize_models)
{
  // we don't want to take ownership of the starting tree, so clone it
  // first and push the clone onto the stack
  pll_utree_t* tree = pll_utree_clone(starting_position.GetTree());
  pll_utree_every(tree, pll::cb_copy_clv_traversal);

  path_.emplace(tree, starting_position.GetModel());
}

Wanderer::~Wanderer()
{
  while (!path_.empty()) {
    pll_utree_destroy(path_.top().GetTree(), pll::cb_erase_data);
    path_.pop();
  }
}

void Wanderer::Start()
{
  //
  // It's important that Start() behaves as if the wanderer arrived at
  // its starting tree during exploration, so we do the same things
  // here that are done in MoveForward(). We check that the tree
  // hasn't been visited, optimize its branch lengths and see if it's
  // good. If this isn't done, teleported wanderers will behave
  // differently than exploring ones.
  //

  pll_utree_t* tree = path_.top().GetTree();

  // set initial model
  partition_.SetModel(path_.top().GetModel());

  if (!authority_.RequestVisit(tree)) {
    // if the starting tree has already been visited, we're done.
    pll_utree_destroy(tree, pll::cb_erase_data);
    path_.pop();

    return;
  }

  // get an inner node of the tree
  pll_unode_t* root = GetVirtualRoot(tree);

  // synchronize partition and tree with a full traversal
  partition_.TraversalUpdate(root, pll::TraversalType::FULL);

  // do optimization. the function will handle its own traversal
  // updates.
  if (optimize_models_) {
    partition_.OptimizeAllBranchesAndModel(root);
  } else {
    partition_.OptimizeAllBranches(root);
  }

  // orient CLVs and compute log-likelihood
  partition_.TraversalUpdate(root, pll::TraversalType::PARTIAL);
  double lnl = partition_.LogLikelihood(root);

  // report the score to the authority. if it returns false, this
  // isn't a good tree and we're done with it, so destroy it and
  // return
  if (!authority_.ReportVisitScore(tree, lnl)) {
    // if the starting tree isn't good, we're done.
    pll_utree_destroy(tree, pll::cb_erase_data);
    path_.pop();

    return;
  }

  QueueMoves();

  // note that move_queues_.top() and move_queues_.top().front() can
  // change during iteration, so we don't store references to them
  // inside the loops. we always want to be looking at the next move
  // in the top queue.
  while (!move_queues_.empty()) {
    while (!move_queues_.top().empty()) {
      MoveForward();
    }

    // we have an empty move queue for the current tree, so move back
    MoveBack();
  }
}

bool Wanderer::TestMove(pll_utree_t* tree, pll_unode_t* node, MoveType type)
{
  bool accept_move;
  double move_score;

  // move testers are in charge of performing whatever traversals are
  // required, and must ensure the partition and tree are synchronized
  // upon return. this does not necessarily mean that all the CLVs are
  // valid and pointed at node, merely that the next partial traversal
  // will work properly.

  std::tie(accept_move, move_score) =
      move_tester_->EvaluateMove(partition_, tree, node, type, authority_);

  // report the test score to the authority
  authority_.ReportTestScore(tree, node, type, move_score);

  return accept_move;
}

void Wanderer::MoveForward()
{
  // the move actually has to be applied to the *original* tree before
  // it's cloned, since the node pointer stored in the move queue
  // points at that tree and not the clone. so we apply the move,
  // clone the tree, and then reverse the move on the original tree.

  // apply move and invalidate the CLVs on that edge
  TreeMove move = move_queues_.top().front();
  move_queues_.top().pop();

  pll_utree_nni(move.node, move.type, nullptr);
  pll::InvalidateEdgeClvs(move.node);

  // clone the tree with move applied (and with properly invalidated CLVs)
  pll_utree_t* tree = pll_utree_clone(path_.top().GetTree());
  pll_utree_every(tree, pll::cb_copy_clv_traversal);

  // undo the move on the original tree. the CLVs on that edge of this tree will
  // remain invalid until the next time it's traversed, which will be either:
  //
  //   * during optimization in the next call to MoveForward(), if
  //     this move doesn't result in a good tree and there are moves
  //     left in the queue; or
  //
  //   * when a call to MoveBack() returns to this tree, at which
  //     point a full traversal is performed if there are moves left
  //     in the queue
  pll_utree_nni(move.node, move.type, nullptr);

  // get an inner node of the cloned tree
  pll_unode_t* root = GetVirtualRoot(tree);

  // do optimization. the function will handle its own traversal
  // updates.
  if (optimize_models_) {
    partition_.OptimizeAllBranchesAndModel(root);
  } else {
    partition_.OptimizeAllBranches(root);
  }

  // orient CLVs and compute log-likelihood
  partition_.TraversalUpdate(root, pll::TraversalType::PARTIAL);
  double lnl = partition_.LogLikelihood(root);

  // report the score to the authority. if it returns true, this is a
  // good tree, and we should push it onto the stack and queue moves
  if (authority_.ReportVisitScore(tree, lnl)) {
    path_.emplace(tree, partition_.GetModel());
    QueueMoves();
  } else {
    // otherwise, we're done with this tree, so destroy it
    pll_utree_destroy(tree, pll::cb_erase_data);
  }
}

void Wanderer::MoveBack()
{
  // this function must perfectly restore the tree and partition to
  // its previous state, so it must ensure that branch lengths are all
  // restored as well as rolling back the NNI move

  // pop the current tree's move queue off the top of the stack
  move_queues_.pop();

  // restore previous position by destroying the current tree and popping
  // off the top of the path stack
  pll_utree_destroy(path_.top().GetTree(), pll::cb_erase_data);
  path_.pop();

  // if the path stack is empty, we're done.
  if (path_.empty()) {
    return;
  }

  // restore the previous model and resynchronize the partition and
  // tree with a full traversal. if the move queue for the tree is
  // empty, this is pointless and can be skipped.
  if (!move_queues_.top().empty()) {
    pll_unode_t* root = GetVirtualRoot(path_.top().GetTree());

    partition_.SetModel(path_.top().GetModel());
    partition_.TraversalUpdate(root, pll::TraversalType::FULL);
  }
}

void Wanderer::QueueMoves()
{
  Position position = path_.top();

  pll_utree_t* tree = position.GetTree();
  pll_unode_t** inner_nodes = tree->nodes + tree->tip_count;

  std::queue<TreeMove> move_queue;
  for (size_t i = 0; i < tree->inner_count; ++i) {
    pll_unode_t* node = inner_nodes[i];

    // skip any pendant edges
    if (!node->back->next) {
      continue;
    }

    // the wanderer doesn't know this, but the way the authority's
    // ProposeMove() method is written, multiple wanderers can be
    // granted permission to test moves to the same tree.
    // ProposeMove() only checks to see if the proposed tree has
    // already been visited, in which case testing the move isn't
    // necessary. this is safe because the authority will still only
    // grant permission to accept the move to a single wanderer, but
    // it's potentially inefficient as moves resulting in the same
    // tree may be tested simultaneously by different wanderers.

    if (authority_.ProposeMove(tree, node, MoveType::LEFT)
        && TestMove(tree, node, MoveType::LEFT)
        && authority_.RequestMove(position, node, MoveType::LEFT)) {
      move_queue.push(TreeMove{node, MoveType::LEFT});
    }

    if (authority_.ProposeMove(tree, node, MoveType::RIGHT)
        && TestMove(tree, node, MoveType::RIGHT)
        && authority_.RequestMove(position, node, MoveType::RIGHT)) {
      move_queue.push(TreeMove{node, MoveType::RIGHT});
    }
  }

  // it's legitimate to push an empty queue onto the stack here; if
  // that happens, it will just result in a move back on the next
  // iteration.
  move_queues_.push(std::move(move_queue));
}

// this function should only be called when path_ and move_queues_
// are empty, i.e., the Wanderer is idle
void Wanderer::Teleport(const Position& position)
{
  if (!path_.empty()) {
    throw std::logic_error("path_ is not empty");
  }

  if (!move_queues_.empty()) {
    throw std::logic_error("move_queues_ is not empty");
  }

  // clone the tree and push it and the model onto the stack
  pll_utree_t* tree = pll_utree_clone(position.GetTree());
  pll_utree_every(tree, pll::cb_copy_clv_traversal);

  path_.emplace(tree, position.GetModel());
}

} // namespace pt
