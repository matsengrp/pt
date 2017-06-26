#ifndef PT_MOVE_TESTER_HPP_
#define PT_MOVE_TESTER_HPP_

#include <utility>

#include <libpll/pll.h>
#include <pll_partition.hpp>

#include "authority.hpp"
#include "wanderer.hpp" // TODO: for MoveType, which should be in another file

namespace pt {

class MoveTester {
 public:
  virtual ~MoveTester();

  virtual std::pair<bool, double> EvaluateMove(pll::Partition& partition,
                                               pll_utree_t* tree,
                                               pll_unode_t* node, MoveType type,
                                               const Authority& authority) const = 0;
};

} // namespace pt

#endif /* PT_MOVE_TESTER_HPP_ */
