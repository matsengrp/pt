#include "pll-utils.hpp"

namespace pt {

bool TreeHealthyAux(pll_utree_t * tree) {
  if (!tree->length)
    return false;

  if (tree->next) {
    if (!tree->next->length || !tree->next->next->length)
      return false;

    return (
        TreeHealthyAux(tree->next->back) &&
          TreeHealthyAux(tree->next->next->back));
  }

  return true;
}


/// @brief Determine if a tree is healthy, i.e. has branch lengths.
/// @param[in] A pll_utree_t.
bool TreeHealthy(pll_utree_t * tree) {
  return (TreeHealthyAux(tree) && TreeHealthyAux(tree->back));
}

}
