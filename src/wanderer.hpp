#ifndef PT_WANDERER_HPP_
#define PT_WANDERER_HPP_

#include <queue>
#include <stack>
#include <string>
#include <vector>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

#include <model_parameters.hpp>
#include <pll_partition.hpp>

#include "authority.hpp"

namespace pt {

enum MoveType : int { LEFT = PLL_UTREE_MOVE_NNI_LEFT,
                      RIGHT = PLL_UTREE_MOVE_NNI_RIGHT };

struct TreeMove {
  pll_utree_t* node;
  MoveType type;
};

class Wanderer {
 private:
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
           pll_utree_t* starting_tree, bool try_all_moves = true);

  Wanderer(Authority& authority,
           pll_utree_t* starting_tree,
           const pll::ModelParameters& model_parameters,
           const std::vector<std::string>& labels,
           const std::vector<std::string>& sequences,
           bool try_all_moves = true);

  ~Wanderer();

  void Teleport(pll_utree_t* starting_tree);
  void Start();

 private:
  void QueueMoves();

  bool TestMove(pll_utree_t* node, MoveType type);

  void MoveForward();
  void MoveBack();
};

} // namespace pt

#endif /* PT_WANDERER_HPP_ */
