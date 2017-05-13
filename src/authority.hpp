#ifndef PT_AUTHORITY_HPP_
#define PT_AUTHORITY_HPP_

#include <atomic>
#include <memory>
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
  // perfect-forwarding factory function; see Effective Modern C++ p. 132
  template <typename... Ts>
  static std::shared_ptr<Authority> Create(Ts&&... params);

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

 private:
  Authority(double ml_lnl, double lnl_offset);
};

template <typename... Ts>
std::shared_ptr<Authority> Authority::Create(Ts&&... params)
{
  // we can't use std::make_shared here, because the constructor is private
  return std::shared_ptr<Authority>(new Authority(std::forward<Ts>(params)...));
}

} // namespace pt

#endif // PT_AUTHORITY_HPP_
