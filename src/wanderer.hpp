#ifndef PT_WANDERER_HPP_
#define PT_WANDERER_HPP_

#include <atomic>
#include <queue>
#include <stack>

#include <cuckoohash_map.hh>
#include <libpll/pll.h>

namespace pt {

using TreeTable = cuckoohash_map<std::string, double>;

enum class MoveType { LEFT = PLL_UTREE_MOVE_NNI_LEFT,
                      RIGHT = PLL_UTREE_MOVE_NNI_RIGHT };

enum class TraversalType { FULL, PARTIAL };

struct TreeMove {
  pll_utree_t* node_;
  MoveType type_;
};

class Authority {
 private:
  std::atomic<double> ml_lnl_;
  const double lambda_;

 public:
  Authority(double ml_lnl, double lambda);

  double GetThreshold() const;

  double GetMaximum() const;
  void SetMaximum(double lnl);
};


class Wanderer {
 private:
  static TreeTable good_trees_;
  static TreeTable all_trees_;

  Authority& authority_;

  pll_partition_t* partition_;
  pll_utree_t* tree_;

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

  std::stack<std::queue<TreeMove>> move_queues_;
  std::stack<pll_utree_t*> trees_;

  std::stack<TreeMove> move_history_;

 public:
  Wanderer(Authority& authority, pll_partition_t* partition,
           pll_utree_t* tree);
  ~Wanderer();

  void Start();

 private:
  void QueueMoves();

  void Move(pll_utree_t* node, MoveType type);
  void MoveBack();

  void Teleport(pll_utree_t* tree);

  void OptimizeBranch(pll_utree_t* node);
  void OptimizeAllBranches();

  double LogLikelihood();

  void Update(pll_utree_t* node, TraversalType type);
};

} // namespace pt

#endif /* PT_WANDERER_HPP_ */
