#ifndef PT_AUTHORITY_HPP_
#define PT_AUTHORITY_HPP_

#include <atomic>
#include <string>

#include <cuckoohash_map.hh>

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
  // TODO: make authority constructor private and add a factory
  //       function that returns a shared pointer
  Authority(double ml_lnl, double lnl_offset);

  // TODO: with multiple threads there can be a race between getting
  //       and setting the maximum, so what we really need is a
  //       SetMaximumIfBetter() function with a proper mutex. see the
  //       comments in Wanderer::MoveForward().
  double GetMaximum() const;
  void SetMaximum(double lnl);

  double GetThreshold() const;

  // TODO: these break encapsulation
  TreeTable& GetVisitedTreeTable();
  TreeTable& GetGoodTreeTable();

  void FilterGoodTreeTable(double lnl_threshold);
  void FilterGoodTreeTable();

  bool InsertVisitedTree(const std::string& newick_str, double lnl);
  bool InsertGoodTree(const std::string& newick_str, double lnl);

  bool UpdateVisitedTree(const std::string& newick_str, double lnl);
  bool UpdateGoodTree(const std::string& newick_str, double lnl);
};

} // namespace pt

#endif // PT_AUTHORITY_HPP_
