#include "authority.hpp"

#include <stdexcept>

namespace pt {

//
// Authority
//

Authority::Authority(double ml_lnl, double lnl_offset) :
    ml_lnl_(ml_lnl),
    lnl_offset_(lnl_offset)
{
  if (lnl_offset_ > 0.0) {
    throw std::invalid_argument("lnl_offset is greater than zero");
  }
}

double Authority::GetMaximum() const
{
  return ml_lnl_;
}

void Authority::SetMaximum(double lnl)
{
  ml_lnl_ = lnl;
}

double Authority::GetThreshold() const
{
  return ml_lnl_ + lnl_offset_;
}

TreeTable& Authority::GetVisitedTreeTable()
{
  return visited_trees_;
}

TreeTable& Authority::GetGoodTreeTable()
{
  return good_trees_;
}

bool Authority::InsertVisitedTree(const std::string& newick_str, double lnl)
{
  return visited_trees_.insert(newick_str, lnl);
}

bool Authority::InsertGoodTree(const std::string& newick_str, double lnl)
{
  return good_trees_.insert(newick_str, lnl);
}

bool Authority::UpdateVisitedTree(const std::string& newick_str, double lnl)
{
  return visited_trees_.update(newick_str, lnl);
}

bool Authority::UpdateGoodTree(const std::string& newick_str, double lnl)
{
  return good_trees_.update(newick_str, lnl);
}

} // namespace pt
