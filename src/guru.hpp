#ifndef PT_GURU_HPP_
#define PT_GURU_HPP_

#include <list>
#include <queue>
#include <string>
#include <vector>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

#include <model_parameters.hpp>
#include <pll_partition.hpp>

#include "authority.hpp"
#include "wanderer.hpp"

namespace pt {

class Guru : public Authority {
 private:
  size_t thread_count_;
  unsigned int tip_node_count_;
  pll::ModelParameters model_parameters_;
  std::vector<std::string> labels_;
  std::vector<std::string> sequences_;
  bool try_all_moves_;

  pll::Partition partition_;
  std::queue<pll_utree_t*> starting_trees_;

  // we use a list here because using a vector would require that
  // wanderers were copy-constructible, which they're not because
  // their partitions can't be copied
  std::list<Wanderer> wanderers_;

 public:
  Guru(double lnl_offset,
       size_t thread_count,
       pll_utree_t* starting_tree,
       unsigned int tip_node_count,
       const pll::ModelParameters& model_parameters,
       const std::vector<std::string>& labels,
       const std::vector<std::string>& sequences,
       bool try_all_moves = true);

  ~Guru() override;

  void AddStartingTree(pll_utree_t* starting_tree);

  void Start();
};

} // namespace pt

#endif /* PT_GURU_HPP_ */
