#ifndef ORDERED_TREE_HPP_INCLUDED
#define ORDERED_TREE_HPP_INCLUDED

#include "ordered-tree.hpp"
#include "partition.hpp"
#include "pll-utils.hpp"

namespace pt {
pll_utree_t *SetNewickRoot(pll_utree_t *tree);
bool SetLabelRoot(std::string label, pll_utree_t *tree, pll_utree_t **root);
void RecursiveOrderedNewick(pll_utree_t *tree);
std::string FindRootNode(pll_utree_t *tree);
pll_utree_t *ToOrderedNewick(pll_utree_t *tree);
std::string ToNewick(pll_utree_t *tree);
std::string ToFullNewick(pll_utree_t *tree);
char *utree_short_newick(pll_utree_t *tree);
static char *newick_utree_recurse(pll_utree_t *tree);
}

#endif // ORDERED-TREE_HPP_INCLUDED
