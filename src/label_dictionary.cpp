#include "label_dictionary.hpp"

#include <limits>
#include <string>
#include <vector>

namespace pt {

LabelDictionary::LabelDictionary()
{ }

LabelDictionary::LabelDictionary(const std::vector<std::string>& labels)
{
  if (labels.empty()) {
    throw std::invalid_argument("labels is empty");
  }

  if (labels.size() >= std::numeric_limits<index_type>::max()) {
    throw std::length_error("index type is not large enough");
  }

  for (size_t i = 0; i < labels.size(); ++i) {
    const auto& label = labels[i];

    if (label_to_index_.find(label) != label_to_index_.end()) {
      throw std::invalid_argument("duplicate label");
    }

    label_to_index_[label] = i;
    index_to_label_[i] = label;
  }
}

std::string LabelDictionary::GetLabel(index_type index) const
{
  auto iter = index_to_label_.find(index);
  if (iter == index_to_label_.end()) {
    throw std::invalid_argument("index not found in dictionary");
  }

  return iter->second;
}

LabelDictionary::index_type
LabelDictionary::GetIndex(const std::string& label) const
{
  auto iter = label_to_index_.find(label);
  if (iter == label_to_index_.end()) {
    throw std::invalid_argument("label not found in dictionary");
  }

  return iter->second;
}

} // namespace pt
