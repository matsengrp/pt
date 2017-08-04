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

pll_unode_t* CreateInnerNode(pll_unode_t* root,
                             const std::string& label)
{
  //
  // done using mallocs so pll_utree_graph_destroy() can be used
  //

  pll_unode_t* node = static_cast<pll_unode_t*>(malloc(sizeof(pll_unode_t)));

  node->back = root;
  node->next = static_cast<pll_unode_t*>(malloc(sizeof(pll_unode_t)));
  node->next->next = static_cast<pll_unode_t*>(malloc(sizeof(pll_unode_t)));
  node->next->next->next = node;

  if (label.empty()) {
    node->label = nullptr;
  } else {
    node->label = static_cast<char*>(malloc(strlen(label.c_str()) + 1));
    strcpy(node->label, label.c_str());
  }

  node->next->label = node->label;
  node->next->next->label = node->label;

  return node;
}

pll_unode_t* CreateTipNode(pll_unode_t* root,
                           const std::string& label)
{
  pll_unode_t* node = static_cast<pll_unode_t*>(malloc(sizeof(pll_unode_t)));

  node->back = root;
  node->next = nullptr;

  if (label.empty()) {
    node->label = nullptr;
  } else {
    node->label = static_cast<char*>(malloc(strlen(label.c_str()) + 1));
    strcpy(node->label, label.c_str());
  }

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
  labels_.clear();

  // true indicates an internal node, false indicates a tip
  bits_.push_back(true);

  if (root->label) {
    labels_.emplace_back(root->label);
  } else {
    labels_.emplace_back("");
  }

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

  if (root->label) {
    labels_.emplace_back(root->label);
  } else {
    labels_.emplace_back("");
  }

  if (root->next) {
    bits_.push_back(true);
    EncodeSubtree(root->next->back);
    EncodeSubtree(root->next->next->back);
  } else {
    bits_.push_back(false);
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
  for (auto label : labels_) {
    ss << delim << "\"" << label << "\"";
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
  std::list<std::string> labels = labels_;

  pll_unode_t* root = Decode(bits, labels);

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
                                    std::list<std::string>& labels)
{
  if (bits.empty()) {
    throw std::invalid_argument("bits is empty");
  }

  bool b = bits.front();
  bits.pop_front();

  if (!b) {
    throw std::invalid_argument("first bit is not set");
  }

  std::string label = labels.front();
  labels.pop_front();

  pll_unode_t* root = CreateInnerNode(nullptr, label);

  root->back = DecodeSubtree(root, bits, labels);
  root->next->back = DecodeSubtree(root->next, bits, labels);
  root->next->next->back = DecodeSubtree(root->next->next, bits, labels);

  return root;
}

pll_unode_t* CompressedTree::DecodeSubtree(pll_unode_t* root,
                                           std::list<bool>& bits,
                                           std::list<std::string>& labels)
{
  if (bits.empty()) {
    // TODO: can this happen? is it a normal part of the algorithm?
    throw std::runtime_error("bits is empty");
  }

  bool b = bits.front();
  bits.pop_front();

  std::string label = labels.front();
  labels.pop_front();

  pll_unode_t* node = nullptr;

  if (b) {
    // inner node
    node = CreateInnerNode(root, label);

    node->next->back = DecodeSubtree(node->next, bits, labels);
    node->next->next->back = DecodeSubtree(node->next->next, bits, labels);
  } else {
    // tip node
    node = CreateTipNode(root, label);
  }

  return node;
}

size_t CompressedTree::Hash() const
{
  size_t seed = std::hash<std::vector<bool>>()(bits_);

  for (auto label : labels_) {
    // from Boost's hash_combine()
    // https://stackoverflow.com/a/20511429
    seed ^= std::hash<std::string>()(label) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  return seed;
}

} // namespace pt
