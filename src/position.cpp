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

pll_utree_t* Position::GetTree() const
{
  return tree_;
}

pll::Model Position::GetModel() const
{
  return model_;
}

} // namespace pt
