#include "compressed_tree.hpp"

#include <cstdlib>
#include <cstring>
#include <functional>
#include <list>
#include <sstream>
#include <vector>

#include <libpll/pll.h>

#include "ordered_tree.hpp"

// TODO: debugging only
#include <iostream>
//#define PRINT_TREES

namespace pt {

//
// free functions
//

pll_unode_t* CreateInnerNode(pll_unode_t* root)
{
  //
  // done using mallocs so pll_utree_graph_destroy() can be used
  //

  pll_unode_t* node = static_cast<pll_unode_t*>(malloc(sizeof(pll_unode_t)));

  node->back = root;
  node->next = static_cast<pll_unode_t*>(malloc(sizeof(pll_unode_t)));
  node->next->next = static_cast<pll_unode_t*>(malloc(sizeof(pll_unode_t)));
  node->next->next->next = node;

  node->label = nullptr;

  node->next->label = node->label;
  node->next->next->label = node->label;

  return node;
}

pll_unode_t* CreateTipNode(pll_unode_t* root,
                           const std::string& label)
{
  if (label.empty()) {
    throw std::invalid_argument("empty tip node label");
  }

  pll_unode_t* node = static_cast<pll_unode_t*>(malloc(sizeof(pll_unode_t)));

  node->back = root;
  node->next = nullptr;

  node->label = static_cast<char*>(malloc(strlen(label.c_str()) + 1));
  strcpy(node->label, label.c_str());

  return node;
}

//
// CompressedTree
//

LabelDictionary CompressedTree::label_dictionary_ = LabelDictionary();

void CompressedTree::BuildDictionary(const std::vector<std::string>& labels)
{
  label_dictionary_ = LabelDictionary(labels);
}

CompressedTree::CompressedTree(const pll_utree_t* tree)
{
  Encode(tree);
}

void CompressedTree::Encode(const pll_utree_t* original_tree)
{
  if (!original_tree) {
    throw std::invalid_argument("tree is null");
  }

  // clone the tree, because ToOrderedNewick() reorders the tree in
  // place and we don't want to modify the tree we were given. note
  // that we aren't modifying any of the node user data (via the
  // node->data pointers) so we don't have to copy or free those.
  pll_utree_t* tree = pll_utree_clone(original_tree);

#ifdef PRINT_TREES
  std::cerr << "original tree\n\n";
  pll_utree_show_ascii(GetVirtualRoot(tree), PLL_UTREE_SHOW_LABEL);
  std::cerr << "\n\n";
#endif

  // ToOrderedNewick() only reorders the tree, despite its name. its
  // return value is the rerooted and reordered tree, which we assign
  // to our pointer before continuing.
  pll_unode_t* root = ToOrderedNewick(GetVirtualRoot(tree));

#ifdef PRINT_TREES
  std::cerr << "reordered tree\n\n";
  pll_utree_show_ascii(root, PLL_UTREE_SHOW_LABEL);
  std::cerr << "\n\n";
#endif

  if (!root) {
    throw std::invalid_argument("root is null");
  }

  if (!root->next) {
    throw std::invalid_argument("Encode() requires an internal node");
  }

  bits_.clear();
  label_indices_.clear();

  // true indicates an internal node, false indicates a tip. we only
  // store labels for tip nodes.
  bits_.push_back(true);

  EncodeSubtree(root->back);
  EncodeSubtree(root->next->back);
  EncodeSubtree(root->next->next->back);

  pll_utree_destroy(tree, nullptr);
}

void CompressedTree::EncodeSubtree(const pll_unode_t* root)
{
  if (!root) {
    throw std::invalid_argument("root is null");
  }

  if (root->next) {
    bits_.push_back(true);
    EncodeSubtree(root->next->back);
    EncodeSubtree(root->next->next->back);
  } else {
    index_type label_index = label_dictionary_.GetIndex(root->label);

    bits_.push_back(false);
    label_indices_.push_back(label_index);
  }
}

std::string CompressedTree::ToDebugString() const
{
  std::stringstream ss;

  ss << "[";
  for (auto bit : bits_) {
    ss << (bit ? "1" : "0");
  }
  ss << "]\n";

  std::string delim = "";

  ss << "[";
  for (auto label_index : label_indices_) {
    ss << delim << label_index;
    delim = ", ";
  }
  ss << "]\n";

  return ss.str();
}

std::string CompressedTree::Decode() const
{
  // the decoding algorithm is destructive, so we make copies of the
  // data members and pass them to static functions that actually do
  // the work. we also can't pop things off the front of a vector, so
  // we copy the bits into a list instead.
  std::list<bool> bits(bits_.begin(), bits_.end());
  std::list<index_type> label_indices = label_indices_;

  pll_unode_t* root = Decode(bits, label_indices);

#ifdef PRINT_TREES
  std::cerr << "decoded tree\n\n";
  pll_utree_show_ascii(root, PLL_UTREE_SHOW_LABEL);
  std::cerr << "\n\n";
#endif

  std::string newick_str = ToNewick(root);

  pll_utree_graph_destroy(root, nullptr);
  return newick_str;
}

pll_unode_t* CompressedTree::Decode(std::list<bool>& bits,
                                    std::list<index_type>& label_indices)
{
  if (bits.empty()) {
    throw std::invalid_argument("bits is empty");
  }

  bool b = bits.front();
  bits.pop_front();

  if (!b) {
    throw std::invalid_argument("first bit is not set");
  }

  pll_unode_t* root = CreateInnerNode(nullptr);

  root->back = DecodeSubtree(root, bits, label_indices);
  root->next->back = DecodeSubtree(root->next, bits, label_indices);
  root->next->next->back = DecodeSubtree(root->next->next, bits, label_indices);

  return root;
}

pll_unode_t* CompressedTree::DecodeSubtree(pll_unode_t* root,
                                           std::list<bool>& bits,
                                           std::list<index_type>& label_indices)
{
  if (bits.empty()) {
    // TODO: can this happen? is it a normal part of the algorithm?
    throw std::runtime_error("bits is empty");
  }

  bool b = bits.front();
  bits.pop_front();

  pll_unode_t* node = nullptr;

  if (b) {
    // inner node
    node = CreateInnerNode(root);

    node->next->back = DecodeSubtree(node->next, bits, label_indices);
    node->next->next->back = DecodeSubtree(node->next->next, bits, label_indices);
  } else {
    // tip node
    std::string label = label_dictionary_.GetLabel(label_indices.front());
    label_indices.pop_front();

    node = CreateTipNode(root, label);
  }

  return node;
}

size_t CompressedTree::Hash() const
{
  size_t seed = std::hash<std::vector<bool>>()(bits_);

  for (auto label_index : label_indices_) {
    // from Boost's hash_combine()
    // https://stackoverflow.com/a/20511429
    seed ^= std::hash<index_type>()(label_index) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  return seed;
}

} // namespace pt
