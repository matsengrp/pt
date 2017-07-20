#ifndef PT_COMMON_HPP_
#define PT_COMMON_HPP_

#include <libpll/pll.h>

namespace pt {

enum MoveType : int { LEFT = PLL_UTREE_MOVE_NNI_LEFT,
      RIGHT = PLL_UTREE_MOVE_NNI_RIGHT };

struct TreeMove {
  pll_unode_t* node;
  MoveType type;
};

} // namespace pt

#endif /* PT_COMMON_HPP_ */
