#ifndef PT_AUTHORITY_HPP_
#define PT_AUTHORITY_HPP_

#include <atomic>
#include <string>
#include <utility>

#include <cuckoohash_map.hh>
#include <libpll/pll.h>

#include "common.hpp"

namespace pt {

using TreeTable = cuckoohash_map<std::string, double>;

//
// Authority
//

class Authority {
 private:
  TreeTable visited_trees_;
  TreeTable good_trees_;

  std::atomic<double> ml_lnl_;
  const double lnl_offset_;

 public:
  Authority(double ml_lnl, double lnl_offset);

  virtual ~Authority();

  double GetMaximumScore() const;
  void SetMaximumScore(double lnl);
  double GetThresholdScore() const;
  std::string GetKey(pll_utree_t* tree) const;

  // TODO: these break encapsulation
  TreeTable& GetVisitedTreeTable();
  TreeTable& GetGoodTreeTable();

  void FilterGoodTreeTable(double lnl_threshold);
  void FilterGoodTreeTable();

  bool ProposeMove(pll_utree_t* tree, pll_unode_t* node, MoveType type);
  virtual bool RequestMove(pll_utree_t* tree, pll_unode_t* node, MoveType type);
  bool RequestTree(pll_utree_t* tree);

  // returns true if the tree is good, false otherwise
  bool ReportTreeScore(pll_utree_t* tree, double lnl);
  bool ReportTreeScore(const std::string& newick_str, double lnl);
};

} // namespace pt

#endif // PT_AUTHORITY_HPP_
