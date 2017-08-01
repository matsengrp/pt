#ifndef PT_MOVE_TESTER_BRANCH_NEIGHBORHOOD_OPTIMIZER_HPP_
#define PT_MOVE_TESTER_BRANCH_NEIGHBORHOOD_OPTIMIZER_HPP_

#include <utility>

#include <libpll/pll.h>
#include <pll_partition.hpp>

#include "../authority.hpp"
#include "../common.hpp"
#include "../move_tester.hpp"

namespace pt { namespace move_tester {

class BranchNeighborhoodOptimizer : public MoveTester
{
 private:
  int radius_;

 public:
  explicit BranchNeighborhoodOptimizer(int radius);
  ~BranchNeighborhoodOptimizer() override;

  std::pair<bool, double> EvaluateMove(pll::Partition& partition,
                                       pll_utree_t* original_tree,
                                       pll_unode_t* original_node,
                                       MoveType type,
                                       const Authority& authority) const override;
};

} } // namespace pt::move_tester

#endif /* PT_MOVE_TESTER_BRANCH_NEIGHBORHOOD_OPTIMIZER_HPP_*/
