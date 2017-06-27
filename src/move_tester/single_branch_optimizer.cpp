#include "single_branch_optimizer.hpp"

#include <utility>

#include <libpll/pll.h>
#include <pll_partition.hpp>

#include "../authority.hpp"
#include "../common.hpp"

namespace pt { namespace move_tester {

SingleBranchOptimizer::~SingleBranchOptimizer()
{ }

std::pair<bool, double>
SingleBranchOptimizer::EvaluateMove(pll::Partition& partition, pll_utree_t* tree,
                                   pll_unode_t* node, MoveType type,
                                   const Authority& authority) const
{
  // apply move and invalidate the CLVs on that edge
  //
  // TODO: should we clone the tree after the move and operate on that
  //       instead, like in Wanderer::MoveForward()?
  pll_utree_nni(node, type, nullptr);
  pll::InvalidateEdgeClvs(node);

  const double original_length = node->length;

  // orient CLVs and optimize branch
  partition.TraversalUpdate(node, pll::TraversalType::PARTIAL);
  partition.OptimizeBranch(node);

  // no need for another traversal, since the CLVs are already
  // pointing at node
  const double test_lnl = partition.LogLikelihood(node);

  // we're just testing whether or not to try the move, so we don't
  // report the score to the authority yet
  bool accept_move = false;
  if (test_lnl >= authority.GetThresholdScore()) {
    accept_move = true;
  }

  // restore the branch length, undo the move, and invalidate the CLVs
  // on that edge again. future operations on this tree will require a
  // traversal first.
  partition.UpdateBranchLength(node, original_length);
  pll_utree_nni(node, type, nullptr);
  pll::InvalidateEdgeClvs(node);

  return std::make_pair(accept_move, test_lnl);
}

} } // namespace pt::move_tester
