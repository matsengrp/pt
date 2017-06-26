#ifndef PT_GURU_HPP_
#define PT_GURU_HPP_

#include <atomic>
#include <deque>
#include <future>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include <libpll/pll.h>
#include <model_parameters.hpp>
#include <pll_partition.hpp>

#include "authority.hpp"
#include "common.hpp"
#include "move_tester.hpp"
#include "wanderer.hpp"

namespace pt {

class Guru : public Authority {
 private:
  size_t thread_count_;
  pll::ModelParameters model_parameters_;
  std::vector<std::string> labels_;
  std::vector<std::string> sequences_;
  std::shared_ptr<const MoveTester> move_tester_;

  pll::Partition partition_;
  pll_utree_t* default_tree_;
  std::queue<pll_utree_t*> starting_trees_;
  std::mutex mutex_;

  // we use a deque here because using a vector would require that
  // wanderers were copy-constructible, which they're not because
  // their partitions can't be copied
  std::deque<Wanderer> wanderers_;
  std::deque<std::future<void>> futures_;
  std::atomic<size_t> idle_wanderer_count_;
  std::vector<bool> wanderer_ready_;

 public:
  Guru(double lnl_offset,
       size_t thread_count,
       pll_utree_t* starting_tree,
       const pll::ModelParameters& model_parameters,
       const std::vector<std::string>& labels,
       const std::vector<std::string>& sequences,
       std::shared_ptr<const MoveTester> move_tester);

  ~Guru() override;

  void AddStartingTree(pll_utree_t* starting_tree);

  bool RequestMove(pll_utree_t* tree, pll_unode_t* node, MoveType type) override;

  void Start();
  void Wait();

 private:
  void UpdateWandererStatus(size_t index);
};

} // namespace pt

#endif /* PT_GURU_HPP_ */
