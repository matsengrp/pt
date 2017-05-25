#include "wanderer.hpp"

#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

#include <model_parameters.hpp>
#include <pll_partition.hpp>
#include <pll_util.hpp>

#include "authority.hpp"

namespace pt {

//
// Wanderer
//

Wanderer::Wanderer(Authority& authority,
                   pll::Partition&& partition, pll_utree_t* starting_tree,
                   bool try_all_moves) :
    authority_(authority),
    partition_(std::move(partition)),
    try_all_moves_(try_all_moves)
{
  // we don't want to take ownership of starting_tree, so clone it
  // first and push the clone onto the stack
  pll_utree_t* tree = pll_utree_clone(starting_tree);
  pll_utree_every(tree, pll::cb_copy_clv_traversal);

  trees_.push(tree);
}

Wanderer::Wanderer(Authority& authority,
                   pll_utree_t* starting_tree, unsigned int tip_node_count,
                   const pll::ModelParameters& model_parameters,
                   const std::vector<std::string>& labels,
                   const std::vector<std::string>& sequences,
                   bool try_all_moves) :
    authority_(authority),
    partition_(starting_tree, tip_node_count, model_parameters, labels, sequences),
    try_all_moves_(try_all_moves)
{
  // we don't want to take ownership of starting_tree, so clone it
  // first and push the clone onto the stack
  pll_utree_t* tree = pll_utree_clone(starting_tree);
  pll_utree_every(tree, pll::cb_copy_clv_traversal);

  // TODO: should we do a TraversalUpdate() here? Start() does one so
  //       we probably don't need to, but it might be good to ensure
  //       that the wanderer is ready as soon as the constructor
  //       returns

  trees_.push(tree);
}

Wanderer::~Wanderer()
{
  while (!trees_.empty()) {
    pll_utree_every(trees_.top(), pll::cb_erase_data);
    pll_utree_destroy(trees_.top());
    trees_.pop();
  }
}

void Wanderer::Start()
{
  // mark the starting tree as visited, then check its log-likelihood
  // (without any branch length optimization) and see if it should be
  // added to the "good" table. if this isn't done, the starting tree
  // may never be visited.
  //
  // TODO: should we optimize the starting tree first?
  pll_utree_t* tree = trees_.top();

  // RequestTree() will also return the ordered Newick string used as
  // the table key that we can use below if the request is accepted.
  bool request_accepted;
  std::string newick_str;
  std::tie(request_accepted, newick_str) = authority_.RequestTree(tree);

  if (!request_accepted) {
    // if the starting tree has already been visited, we're done.
    return;
  }

  partition_.TraversalUpdate(tree, pll::TraversalType::FULL);
  double lnl = partition_.LogLikelihood(tree);

  // report the score to the authority. we don't care if it's good or
  // not, as we're going to queue the available moves anyway.
  authority_.ReportTreeScore(newick_str, lnl);

  QueueMoves();

  // note that move_queues_.top() and move_queues_.top().front() can
  // change during iteration, so we don't store references to them
  // inside the loops. we always want to be looking at the next move
  // in the top queue.
  while (!move_queues_.empty()) {
    while (!move_queues_.top().empty()) {
      // somewhere here is where we would check to see if there are
      // any idle monks, and use Teleport() (or ask the authority to)
      // to have the idle monk go there instead of taking the move
      // ourselves

      MoveForward();
    }

    // we have an empty move queue for the current tree, so move back
    MoveBack();
  }
}

bool Wanderer::TestMove(pll_utree_t* node, MoveType type)
{
  pll_utree_nni(node, type, nullptr);

  // request permission to proceed from the authority. if successful,
  // this Wanderer owns the tree and we can proceed. RequestTree()
  // will also return the ordered Newick string used as the table key
  // that we can use below if the request is accepted.
  bool request_accepted;
  std::string newick_str;
  std::tie(request_accepted, newick_str) = authority_.RequestTree(node);

  if (!request_accepted) {
    // undo the move and reject
    pll_utree_nni(node, type, nullptr);
    return false;
  }

  //
  // If try_all_moves_ is true, TestMove() will always return true if
  // this tree hasn't been visited before. Otherwise, TestMove() will
  // optimize the branch length and use the resulting log-likelihood
  // to determine if this move should be accepted.
  //

  if (try_all_moves_) {
    // undo the move and accept
    pll_utree_nni(node, type, nullptr);
    return true;
  } else {
    const double original_length = node->length;

    // update so that CLVs are pointing at node and optimize branch
    partition_.TraversalUpdate(node, pll::TraversalType::PARTIAL);
    partition_.OptimizeBranch(node);

    // when OptimizeBranch() returns, where will the CLVs be oriented?
    // do CLVs need to be invalidated? my belief is that if we use the
    // optimized edge for computing the log-likelihood, no CLVs need to
    // be invalidated first.

    // no need for another traversal, since the CLVs are already
    // pointing at node
    const double test_lnl = partition_.LogLikelihood(node);

    // we're just testing whether or not to try the move, so we don't
    // report the score to the authority yet
    bool accept_move = false;
    if (test_lnl >= authority_.GetThresholdScore()) {
      accept_move = true;
    }

    // restore the branch length and undo the move. CLVs will remain
    // pointed toward node, so future operations on this tree will need
    // to orient the CLVs as appropriate
    partition_.UpdateBranchLength(node, original_length);
    pll_utree_nni(node, type, nullptr);

    return accept_move;
  }
}

void Wanderer::MoveForward()
{
  // the move actually has to be applied to the *original* tree before
  // it's cloned, since the node pointer stored in the move queue
  // points at that tree and not the clone. so we apply the move,
  // clone the tree, and then reverse the move on the original tree.

  // apply move. the move we're applying should be
  // move_queues_.top().front(), I think
  TreeMove move = move_queues_.top().front();
  move_queues_.top().pop();

  pll_utree_nni(move.node, move.type, nullptr);

  // clone the tree with move applied
  pll_utree_t* tree = pll_utree_clone(trees_.top());
  pll_utree_every(tree, pll::cb_copy_clv_traversal);

  // undo the move on the original tree
  pll_utree_nni(move.node, move.type, nullptr);

  // do full branch optimization. this function will handle its own
  // traversal updates.
  partition_.OptimizeAllBranches(tree);

  // orient CLVs and compute log-likelihood
  partition_.TraversalUpdate(tree, pll::TraversalType::PARTIAL);
  double lnl = partition_.LogLikelihood(tree);

  // report the score to the authority. if it returns true, this is a
  // good tree, and we should push it onto the stack and queue moves
  if (authority_.ReportTreeScore(tree, lnl)) {
    trees_.push(tree);
    QueueMoves();
  } else {
    // otherwise, we're done with this tree, so destroy it
    pll_utree_every(tree, pll::cb_erase_data);
    pll_utree_destroy(tree);
  }
}

void Wanderer::MoveBack()
{
  // this function must perfectly restore the tree and partition to
  // its previous state, so it must ensure that branch lengths are all
  // restored as well as rolling back the NNI move

  // pop the current tree's move queue off the top of the stack
  move_queues_.pop();

  // restore previous tree by destroying the current tree and popping
  // it off the top of the tree history stack
  pll_utree_every(trees_.top(), pll::cb_erase_data);
  pll_utree_destroy(trees_.top());
  trees_.pop();

  // if the tree stack is empty, we're done.
  if (trees_.empty()) {
    return;
  }

  // resynchronize partition and tree with a full traversal. if the
  // move queue for the tree is empty, this is pointless and can be
  // skipped.
  if (!move_queues_.top().empty()) {
    partition_.TraversalUpdate(trees_.top(), pll::TraversalType::FULL);
  }
}

void Wanderer::QueueMoves()
{
  // is n_tips - 2 correct for the number of inner nodes? probably
  // should add an error check to see if the number of nodes
  // pll_utree_query_innernodes() finds is the same as the size of the
  // vector
  std::vector<pll_utree_t*> inner_nodes(partition_.inner_node_count(), nullptr);
  pll_utree_query_innernodes(trees_.top(), inner_nodes.data());

  std::queue<TreeMove> move_queue;
  for (auto node : inner_nodes) {
    // skip any pendant edges
    if (!node->back->next) {
      continue;
    }

    // TODO: I think that the traversal that starts at trees_.top()
    //       produces a collection of inner nodes that includes the
    //       node at the other end of the "current" edge (that is,
    //       trees_.top()->back). this results in testing duplicate
    //       moves across that edge, doesn't it? we might need
    //       something like the following to avoid this
    //
    // if (node == trees_.top()->back) {
    //   continue;
    // }

    // here is where we can try applying a move, optimizing the
    // current branch length, and then comparing the log-likelihood to
    // the threshold to see if the move is worth taking (since if we
    // do take the move, we'll be performing a full branch length
    // optimization). if the move is worth taking, we add the move to
    // the move queue. regardless, after testing the move, we undo the
    // move and restore the original branch length, then repeat for
    // the move in the other "direction".

    // assume we have a TestMove() method that applies the move,
    // optimizes the branch length, and compares the log-likelihood to
    // the authority's threshold. it then reverses the move and
    // restores the branch length. it may be okay if the CLVs are
    // oriented differently than they were before TestMove() gets
    // called, since TestMove() will have to do a partial traversal
    // anyway, right? TestMove() will return true if the move is worth
    // taking, false otherwise.

    if (TestMove(node, MoveType::LEFT)) {
      move_queue.push(TreeMove{node, MoveType::LEFT});
    }

    if (TestMove(node, MoveType::RIGHT)) {
      move_queue.push(TreeMove{node, MoveType::RIGHT});
    }
  }

  // it's legitimate to push an empty queue onto the stack here; if
  // that happens, it will just result in a move back on the next
  // iteration.
  move_queues_.push(std::move(move_queue));
}

// this function should only be called when trees_ and move_queues_
// are empty, i.e., the Wanderer is idle
void Wanderer::Teleport(pll_utree_t* starting_tree)
{
  if (!trees_.empty()) {
    throw std::logic_error("trees_ is not empty");
  }

  if (!move_queues_.empty()) {
    throw std::logic_error("move_queues_ is not empty");
  }

  // clone the tree and push it onto the stack
  pll_utree_t* tree = pll_utree_clone(starting_tree);
  pll_utree_every(tree, pll::cb_copy_clv_traversal);

  trees_.push(tree);
}

} // namespace pt
