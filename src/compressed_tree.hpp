#ifndef PT_COMPRESSED_TREE_HPP_
#define PT_COMPRESSED_TREE_HPP_

#include <list>
#include <string>
#include <vector>

#include <libpll/pll.h>

#include "label_dictionary.hpp"

namespace pt {

class CompressedTree {
 private:
  static LabelDictionary label_dictionary_;

  std::vector<bool> bits_;
  std::list<std::string> labels_;

 public:
  static void BuildDictionary(const std::vector<std::string>& labels);

  CompressedTree(const pll_utree_t* original_tree);

  std::string Decode() const;

  bool operator==(const CompressedTree& rhs) const;
  bool operator!=(const CompressedTree& rhs) const;

  std::string ToDebugString() const;

  size_t Hash() const;

 private:
  void Encode(const pll_utree_t* original_tree);
  void EncodeSubtree(const pll_unode_t* root);

  static pll_unode_t* Decode(std::list<bool>& bits,
                             std::list<std::string>& labels);
  static pll_unode_t* DecodeSubtree(pll_unode_t* root,
                                    std::list<bool>& bits,
                                    std::list<std::string>& labels);
};

inline bool CompressedTree::operator==(const CompressedTree& rhs) const
{
  return bits_ == rhs.bits_ && labels_ == rhs.labels_;
}

inline bool CompressedTree::operator!=(const CompressedTree& rhs) const
{
  return !(*this == rhs);
}

} // namespace pt


#endif /* PT_COMPRESSED_TREE_HPP_ */
