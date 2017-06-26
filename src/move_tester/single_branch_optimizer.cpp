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
  // apply the move
  pll_utree_nni(node, type, nullptr);

  const double original_length = node->length;

  // the NNI move invalidated the CLVs at either end of this edge,
  // but the CLV validity flags associated with those nodes are
  // unchanged. perform a full traversal to recompute those CLVs and
  // reach a predictable state.
  //
  // TODO: a full traversal is overkill -- perhaps we could use an
  //       InvalidateEdge() method so that we could still use a
  //       partial traversal?
  partition.TraversalUpdate(node, pll::TraversalType::FULL);
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

  // restore the branch length and undo the move. CLVs will remain
  // pointed toward node, so future operations on this tree will need
  // to orient the CLVs as appropriate
  partition.UpdateBranchLength(node, original_length);
  pll_utree_nni(node, type, nullptr);

  // TODO: note that undoing the NNI move again invalidates the CLVs
  //       at either end of this edge, so future partial traversals
  //       will give a bad result. this is another place where an
  //       InvalidateEdge() method would be helpful. an alternative
  //       would be to clone the input tree before applying the NNI
  //       move and testing the single-branch optimization, so that
  //       we wouldn't have to restore the input tree's state.

  return std::make_pair(accept_move, test_lnl);
}

} } // namespace pt::move_tester
