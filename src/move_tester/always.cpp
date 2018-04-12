#include "always.hpp"

#include <utility>

#include <libpll/pll.h>
#include <pll_partition.hpp>

#include "../authority.hpp"
#include "../common.hpp"
#include "../options.hpp"

namespace pt { namespace move_tester {

Always::Always(const Options& options) :
    MoveTester(options)
{ }

Always::~Always()
{ }

std::pair<bool, double>
Always::EvaluateMove(pll::Partition&, pll_utree_t*,
                     pll_unode_t*, MoveType, const Authority&) const
{
  return std::make_pair(true, 0.0);
}

} } // namespace pt::move_tester
