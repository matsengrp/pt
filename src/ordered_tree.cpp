#include "ordered_tree.hpp"

#include <string>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

namespace pt {

pll_unode_t* GetVirtualRoot(pll_utree_t* tree)
{
  return tree->nodes[tree->tip_count + tree->inner_count - 1];
}

/// @brief Recursive function to generate newick string without branch lengths.
/// @return Newick char* string.
char *newick_utree_recurse(pll_unode_t *tree) {
  char *newick;
  int size_alloced;
  assert(tree != NULL);
  if (!tree->next)
    size_alloced = asprintf(&newick, "%s", tree->label);
  else {
    char *subtree1 = newick_utree_recurse(tree->next->back);
    if (subtree1 == NULL)
      return NULL;
    char *subtree2 = newick_utree_recurse(tree->next->next->back);
    if (subtree2 == NULL) {
      free(subtree1);
      return NULL;
    }
    size_alloced = asprintf(&newick, "(%s,%s)%s", subtree1, subtree2,
                            tree->label ? tree->label : "");
    free(subtree1);
    free(subtree2);
  }
  if (size_alloced < 0) {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
    return NULL;
  }

  return newick;
}

/// @brief Build newick string without branch lengths.
/// @return Newick char* string without branch lengths.
char *utree_short_newick(pll_unode_t *tree) {
  char *newick;
  int size_alloced;

  if (!tree)
    return NULL;
  // If we are given a leaf then move to the attached internal node.
  if (!tree->next)
    tree = tree->back;

  char *subtree1 = newick_utree_recurse(tree->back);
  if (subtree1 == NULL)
    return NULL;
  char *subtree2 = newick_utree_recurse(tree->next->back);
  if (subtree2 == NULL) {
    free(subtree1);
    return NULL;
  }
  char *subtree3 = newick_utree_recurse(tree->next->next->back);
  if (subtree3 == NULL) {
    free(subtree1);
    free(subtree2);
    return NULL;
  }

  size_alloced = asprintf(&newick, "(%s,%s,%s)%s;", subtree1, subtree2,
                          subtree3, tree->label ? tree->label : "");
  free(subtree1);
  free(subtree2);
  free(subtree3);
  if (size_alloced < 0) {
    pll_errno = PLL_ERROR_MEM_ALLOC;
    snprintf(pll_errmsg, 200, "memory allocation during newick export failed");
    return NULL;
  }
  return (newick);
}

/// @brief Returns a the tree as a Newick string.
/// @return Newick string.
std::string ToNewick(pll_unode_t *tree) {
  char *newick = utree_short_newick(tree);

  std::string strnewick;
  if (newick) {
    strnewick = newick;
    free(newick);
  }

  return strnewick;
}

/// @brief Returns a the tree as a Newick string with branch lengths.
/// @return Newick std::string.
std::string ToFullNewick(pll_unode_t *tree) {
  char *newick = pll_utree_export_newick(tree, nullptr);

  std::string strnewick;
  if (newick) {
    strnewick = newick;
    free(newick);
  }

  return strnewick;
}

/// @brief Finds node tree with the label that is first
/// alphabetically, and reroots the tree at its parent node.
/// @param[in,out] tree
/// Tree.
/// @return Smallest label.
std::string FindRootNode(pll_unode_t *tree) {
  std::string minlabel;
  if (!tree->next)
    return tree->label;
  std::string rlabel1 = FindRootNode(tree->next->back);
  std::string rlabel2 = FindRootNode(tree->next->next->back);
  if (rlabel1.compare(rlabel2) < 0)
    minlabel = rlabel1;
  else {
    minlabel = rlabel2;
    pll_unode_t *temp = tree->next;
    tree->next = tree->next->next;
    tree->next->next = temp;
    tree->next->next->next = tree;
  }

  return minlabel;
}

/// @brief Finds a node in a tree with the given label and sets the current node
/// to the
/// parent of it.
/// @param[in] label
/// Label of node to find.
/// @param[in] tree
/// Tree.
/// @param[out] root
/// Buffer for storing root node.
/// @return Were we successful in finding the label?
bool SetLabelRoot(std::string label, pll_unode_t *tree, pll_unode_t **root) {
  if (!tree->next) {
    if (tree->label == label) {
      *root = tree->back;
      return true;
    } else
      return false;
  } else {
    return (SetLabelRoot(label, tree->next->back, root) ||
            SetLabelRoot(label, tree->next->next->back, root));
  }
}

/// @brief Returns the passed tree rooted at the minimum label.
/// @param[in] tree
/// Tree.
/// @return Parent of the alphabetically smallest leaf.
pll_unode_t *SetNewickRoot(pll_unode_t *tree) {
  std::string minlabel;
  std::string rlabel1 = FindRootNode(tree);
  std::string rlabel2 = FindRootNode(tree->back);
  if (rlabel1.compare(rlabel2) < 0)
    minlabel = rlabel1;
  else
    minlabel = rlabel2;
  pll_unode_t *root;
  SetLabelRoot(minlabel, tree, &root);
  SetLabelRoot(minlabel, tree->back, &root);
  return root;
}

/// @brief Determines order for first ternary step, then runs recursive ordering
/// for the two non-root subtrees.
/// @return Completely ordered newick string.
/// @todo Update.
pll_unode_t *ToOrderedNewick(pll_unode_t *tree) {
  tree = SetNewickRoot(tree);
  FindRootNode(tree);
  return tree;
}
}
