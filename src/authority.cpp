#include "authority.hpp"

#include <stdexcept>
#include <string>
#include <vector>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

#include "ordered_tree.hpp"

namespace pt {

//
// free functions
//

std::string OrderedNewickString(pll_utree_t* tree)
{
  // clone the tree, because ToOrderedNewick() reorders the tree in
  // place and we don't want to modify the tree we were given. note
  // that we aren't modifying any of the node user data (via the
  // node->data pointers) so we don't have to copy or free those.
  pll_utree_t* clone = pll_utree_clone(tree);

  // ToOrderedNewick() only reorders the tree, despite its name. its
  // return value is the rerooted and reordered tree, which we assign
  // to our pointer before continuing.
  clone = ToOrderedNewick(clone);
  std::string newick_str = ToNewick(clone);

  pll_utree_destroy(clone);
  return newick_str;
}

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

Authority::~Authority()
{ }

double Authority::GetMaximumScore() const
{
  return ml_lnl_;
}

void Authority::SetMaximumScore(double lnl)
{
  ml_lnl_ = lnl;
}

double Authority::GetThresholdScore() const
{
  return ml_lnl_ + lnl_offset_;
}

std::string Authority::GetKey(pll_utree_t* tree) const
{
  return OrderedNewickString(tree);
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
  FilterGoodTreeTable(GetThresholdScore());
}

std::pair<bool, std::string> Authority::RequestTree(pll_utree_t* tree,
                                                    bool /* first_tree */)
{
  std::string newick_str = GetKey(tree);

  return std::make_pair(visited_trees_.insert(newick_str, 0.0), newick_str);
}

bool Authority::ReportTreeScore(pll_utree_t* tree, double lnl)
{
  return ReportTreeScore(GetKey(tree), lnl);
}

bool Authority::ReportTreeScore(const std::string& newick_str, double lnl)
{
  // TODO: add locks here where appropriate, such as between getting
  //       and setting the maximum

  if (!visited_trees_.contains(newick_str)) {
    throw std::logic_error("tree is not in the visited table");
  }

  if (lnl < GetThresholdScore()) {
    return false;
  }

  if (!good_trees_.insert(newick_str, lnl)) {
    throw std::logic_error("tree is already in the good table");
  }

  if (lnl > ml_lnl_) {
    ml_lnl_ = lnl;
  }

  return true;
}

} // namespace pt
