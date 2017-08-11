#ifndef PT_LABEL_DICTIONARY_HPP_
#define PT_LABEL_DICTIONARY_HPP_

#include <cstdint>
#include <map>
#include <string>
#include <vector>

namespace pt {

class LabelDictionary {
 public:
  using index_type = uint8_t;

 private:
  std::map<std::string, index_type> label_to_index_;
  std::map<index_type, std::string> index_to_label_;

 public:
  LabelDictionary();
  explicit LabelDictionary(const std::vector<std::string>& labels);

  std::string GetLabel(index_type index) const;
  index_type GetIndex(const std::string& label) const;
};

} // namespace pt

#endif /* PT_LABEL_DICTIONARY_HPP_ */
