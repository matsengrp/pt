#ifndef PT_POSITION_HPP_
#define PT_POSITION_HPP_

#include <libpll/pll.h>
#include <model.hpp>

namespace pt {

class Position {
  // TODO: the position class is a very thin wrapper around a tree and
  //       a model and relies on the user for memory management.

 private:
  pll_utree_t* tree_;
  pll::Model model_;

 public:
  Position();
  Position(pll_utree_t* tree, const pll::Model& model);

  pll_utree_t* GetTree() const;
  pll::Model GetModel() const;
};

} // namespace pt

#endif /* PT_POSITION_HPP_ */
