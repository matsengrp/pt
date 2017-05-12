#include "authority.hpp"

#include <stdexcept>
#include <string>
#include <vector>

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

void Authority::FilterGoodTreeTable(double lnl_threshold)
{
  // from the documentation it is unclear if erasing an item from a
  // libcuckoo map invalidates other iterators, so to be safe, we do
  // this in two passes.

  std::vector<std::string> keys_to_erase;

  // lock scope
  {
    auto locked_table = good_trees_.lock_table();
    for (auto& item : locked_table) {
      if (item.second < lnl_threshold) {
        keys_to_erase.push_back(item.first);
      }
    }
  }

  // TODO: in recent versions of libcuckoo, locked tables have their
  //       own erase function. using that, we could guarantee that the
  //       table is filtered properly the moment this function
  //       returns; as is, it's possible for other threads to be
  //       inserting additional items while we're erasing.
  for (auto& key : keys_to_erase) {
    good_trees_.erase(key);
  }
}

void Authority::FilterGoodTreeTable()
{
  FilterGoodTreeTable(GetThreshold());
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
