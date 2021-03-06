#ifndef PT_COMPRESSED_TREE_HPP_
#define PT_COMPRESSED_TREE_HPP_

#include <iostream>
#include <list>
#include <string>
#include <vector>

#include <libpll/pll.h>

#include "label_dictionary.hpp"

namespace pt {

class CompressedTree {
 private:
  using index_type = LabelDictionary::index_type;

  static LabelDictionary label_dictionary_;

  std::vector<bool> bits_;
  std::vector<index_type> label_indices_;

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
                             std::list<index_type>& label_indices);
  static pll_unode_t* DecodeSubtree(pll_unode_t* root,
                                    std::list<bool>& bits,
                                    std::list<index_type>& label_indices);
};

inline bool CompressedTree::operator==(const CompressedTree& rhs) const
{
  return bits_ == rhs.bits_ && label_indices_ == rhs.label_indices_;
}

inline bool CompressedTree::operator!=(const CompressedTree& rhs) const
{
  return !(*this == rhs);
}

inline std::ostream& operator<<(std::ostream& os, const CompressedTree& ct)
{
  os << ct.Decode();
  return os;
}

} // namespace pt

namespace std {

// specialization of std::hash for compressed trees
template <>
struct hash<pt::CompressedTree> {
  inline size_t operator()(const pt::CompressedTree& ct) const
  {
    return ct.Hash();
  }
};

} // namespace std

#endif /* PT_COMPRESSED_TREE_HPP_ */
