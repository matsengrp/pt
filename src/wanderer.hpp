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

namespace pt {

using TreeTable = cuckoohash_map<std::string, double>;

enum class MoveType { LEFT = PLL_UTREE_MOVE_NNI_LEFT,
                      RIGHT = PLL_UTREE_MOVE_NNI_RIGHT };

enum class TraversalType { FULL, PARTIAL };

struct TreeMove {
  pll_utree_t* node;
  MoveType type;
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
  Wanderer(Authority& authority, pll_partition_t* partition,
           pll_utree_t* initial_tree, bool try_all_moves = false);
  ~Wanderer();

  void Start();

 private:
  void QueueMoves();

  bool TestMove(pll_utree_t* node, MoveType type);

  void MoveForward();
  void MoveBack();

  void Teleport(pll_utree_t* tree);

  void OptimizeBranch(pll_utree_t* node);
  void OptimizeAllBranchesOnce(pll_utree_t* tree);
  void OptimizeAllBranches(pll_utree_t* tree);

  double LogLikelihood(pll_utree_t* tree);

  void TraversalUpdate(pll_utree_t* root, TraversalType type);
};

} // namespace pt

#endif /* PT_WANDERER_HPP_ */
