#ifndef PT_WANDERER_HPP_
#define PT_WANDERER_HPP_

#include <atomic>
#include <queue>
#include <stack>

#include <cuckoohash_map.hh>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

#include "pll_partition.hpp"
#include "pll-utils.hpp"

namespace pt {

using TreeTable = cuckoohash_map<std::string, double>;

enum MoveType : int { LEFT = PLL_UTREE_MOVE_NNI_LEFT,
                      RIGHT = PLL_UTREE_MOVE_NNI_RIGHT };

struct TreeMove {
  pll_utree_t* node;
  MoveType type;
};

class Authority {
 private:
  std::atomic<double> ml_lnl_;
  const double lambda_;

 public:
  // TODO: make authority constructor private and add a factory
  //       function that returns a shared pointer
  Authority(double ml_lnl, double lambda);

  double GetThreshold() const;

  // TODO: with multiple threads there can be a race between getting
  //       and setting the maximum, so what we really need is a
  //       SetMaximumIfBetter() function with a proper mutex. see the
  //       comments in Wanderer::MoveForward().
  double GetMaximum() const;
  void SetMaximum(double lnl);

  // TODO: move tree tables here

  // TODO: for when the tree tables are moved here:
  // bool InsertVisitedTree(const std::string& newick_str);
  // bool InsertGoodTree(const std::string& newick_str);
};

class Wanderer {
 private:
  static TreeTable all_trees_;
  static TreeTable good_trees_;

  // TODO: the authority should be a shared pointer constructed via a
  //       factory function
  Authority& authority_;
  pll::Partition partition_;

  const bool try_all_moves_;

  // a data structure that makes more sense for this might be a stack
  // of pairs like {tree, move_queue}. MoveBack() could pop the top of
  // the stack and reset things for the new tree on the top of the
  // stack and continue with the moves in the queue. life would be
  // easier if there were a Tree class that took care of tree resource
  // management, copying etc.

  // the primary benefit of keeping the tree history around is so we
  // don't have to worry about restoring branch lengths upon moving
  // back after a full optimization is performed. I don't think
  // storing all the branch lengths would use much less memory than
  // storing the whole tree.

  std::stack<pll_utree_t*> trees_;
  std::stack<std::queue<TreeMove>> move_queues_;

 public:
  Wanderer(Authority& authority, pll::Partition&& partition,
           pll_utree_t* initial_tree, bool try_all_moves = false);
  ~Wanderer();

  void Start();

  TreeTable& GetVisitedTreeTable();
  TreeTable& GetGoodTreeTable();

 private:
  void QueueMoves();

  bool TestMove(pll_utree_t* node, MoveType type);

  void MoveForward();
  void MoveBack();

  void Teleport(pll_utree_t* tree);
};

inline TreeTable& Wanderer::GetVisitedTreeTable()
{
  return all_trees_;
}

inline TreeTable& Wanderer::GetGoodTreeTable()
{
  return good_trees_;
}

} // namespace pt

#endif /* PT_WANDERER_HPP_ */
