#ifndef PT_MOVE_TESTER_SINGLE_BRANCH_OPTIMIZER_
#define PT_MOVE_TESTER_SINGLE_BRANCH_OPTIMIZER_

#include <utility>

#include <libpll/pll.h>
#include <pll_partition.hpp>

#include "../authority.hpp"
#include "../common.hpp"
#include "../move_tester.hpp"

namespace pt { namespace move_tester {

class SingleBranchOptimizer : public MoveTester
{
 public:
  ~SingleBranchOptimizer() override;

  std::pair<bool, double> EvaluateMove(pll::Partition& partition,
                                       pll_utree_t* tree,
                                       pll_unode_t* node, MoveType type,
                                       const Authority& authority) const override;
};

} } // namespace pt::move_tester

#endif /* PT_MOVE_TESTER_SINGLE_BRANCH_OPTIMIZER_ */
