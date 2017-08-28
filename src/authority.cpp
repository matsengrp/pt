#include "authority.hpp"

#include <stdexcept>
#include <string>
#include <vector>

#include <libpll/pll.h>

#include "common.hpp"
#include "compressed_tree.hpp"
#include "options.hpp"
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
  pll_unode_t* root = ToOrderedNewick(GetVirtualRoot(clone));
  std::string newick_str = ToNewick(root);

  pll_utree_destroy(clone, nullptr);
  return newick_str;
}

//
// Authority
//

Authority::Authority(const Options& options, double ml_lnl) :
    ml_lnl_(ml_lnl),
    lnl_offset_(options.lnl_offset),
    track_tested_trees_(options.track_tested_trees)
{
  if (lnl_offset_ >= 0.0) {
    throw std::invalid_argument("lnl_offset must be less than 0");
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

CompressedTree Authority::GetKey(const pll_utree_t* tree) const
{
  return CompressedTree(tree);
}

TreeTable& Authority::GetTestedTreeTable()
{
  return tested_trees_;
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

  std::vector<CompressedTree> keys_to_erase;

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

bool Authority::ProposeMove(pll_utree_t* tree, pll_unode_t* node, MoveType type)
{
  // apply the move
  pll_utree_nni(node, type, nullptr);

  // check to see if the proposed tree has already been visited
  CompressedTree key = GetKey(tree);
  bool proposal_accepted = !visited_trees_.contains(key);

  // undo the move
  pll_utree_nni(node, type, nullptr);

  return proposal_accepted;
}

bool Authority::RequestMove(pll_utree_t* tree, pll_unode_t* node, MoveType type)
{
  // apply the move
  pll_utree_nni(node, type, nullptr);

  // check to see if the requested tree has already been visited
  CompressedTree key = GetKey(tree);
  bool request_accepted = visited_trees_.insert(key, 0.0);

  // undo the move
  pll_utree_nni(node, type, nullptr);

  return request_accepted;
}

bool Authority::RequestTree(pll_utree_t* tree)
{
  CompressedTree key = GetKey(tree);
  return visited_trees_.insert(key, 0.0);
}

void Authority::ReportTestScore(pll_utree_t* tree, pll_unode_t* node,
                                MoveType type, double score)
{
  if (!track_tested_trees_) {
    return;
  }

  // apply the move
  pll_utree_nni(node, type, nullptr);

  // try and insert the tree into the tested trees table. if it's
  // already there, update it if this score is better.

  // the lambda returns false to indicate this element should not be erased
  auto updater = [score](double& value) {
    if (score > value) { value = score; }
    return false;
  };

  CompressedTree key = GetKey(tree);
  tested_trees_.upsert(key, updater, score);

  // undo the move
  pll_utree_nni(node, type, nullptr);
}

bool Authority::ReportVisitScore(pll_utree_t* tree, double lnl)
{
  return ReportVisitScore(GetKey(tree), lnl);
}

bool Authority::ReportVisitScore(const CompressedTree& key, double lnl)
{
  if (!visited_trees_.update(key, lnl)) {
    throw std::logic_error("tree is not in the visited table");
  }

  if (lnl < GetThresholdScore()) {
    return false;
  }

  if (!good_trees_.insert(key, lnl)) {
    throw std::logic_error("tree is already in the good table");
  }

  {
    std::lock_guard<std::mutex> lock(mutex_);

    if (lnl > ml_lnl_) {
      SetMaximumScore(lnl);
    }
  }

  return true;
}

} // namespace pt
