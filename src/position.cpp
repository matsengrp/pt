#include "position.hpp"

#include <libpll/pll.h>
#include <model.hpp>
#include <pll_util.hpp>

namespace pt {

Position::Position() :
    tree_(nullptr)
{ }

Position::Position(pll_utree_t* tree, const pll::Model& model) :
    tree_(tree),
    model_(model)
{ }

} // namespace pt
