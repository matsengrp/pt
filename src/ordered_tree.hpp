#ifndef ORDERED_TREE_HPP_INCLUDED
#define ORDERED_TREE_HPP_INCLUDED

#include <string>

#include <libpll/pll.h>

namespace pt {

pll_unode_t* GetVirtualRoot(pll_utree_t* tree);

pll_unode_t *SetNewickRoot(pll_unode_t *tree);
bool SetLabelRoot(std::string label, pll_unode_t *tree, pll_unode_t **root);
void RecursiveOrderedNewick(pll_unode_t *tree);
std::string FindRootNode(pll_unode_t *tree);
pll_unode_t *ToOrderedNewick(pll_unode_t *tree);
std::string ToNewick(pll_unode_t *tree);
std::string ToFullNewick(pll_unode_t *tree);
char *utree_short_newick(pll_unode_t *tree);
static char *newick_utree_recurse(pll_unode_t *tree);
}

#endif // ORDERED_TREE_HPP_INCLUDED
