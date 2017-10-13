#include "position.hpp"

#include <libpll/pll.h>
#include <model.hpp>
#include <pll_util.hpp>

namespace pt {

Position::Position(pll_utree_t* tree, pll::Model model) :
    tree_(tree),
    model_(model)
{ }

Position::~Position()
{
  pll_utree_destroy(tree_, pll::cb_erase_data);
}

} // namespace pt
