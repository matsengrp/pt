#include "always.hpp"

#include <utility>

#include <libpll/pll.h>
#include <pll_partition.hpp>

#include "../wanderer.hpp" // TODO: for MoveType, which should be in another file

namespace pt { namespace move_tester {

Always::~Always()
{ }

std::pair<bool, double>
Always::EvaluateMove(pll::Partition&, pll_utree_t*,
                     pll_unode_t*, MoveType, double) const
{
  return std::make_pair(true, 0.0);
}

} } // namespace pt::move_tester
