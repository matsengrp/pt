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
#include <model.hpp>
#include <pll_partition.hpp>

#include "authority.hpp"
#include "common.hpp"
#include "move_tester.hpp"
#include "options.hpp"
#include "position.hpp"
#include "wanderer.hpp"

namespace pt {

class Guru : public Authority {
 private:
  size_t thread_count_;
  pll::Model model_;
  std::vector<std::string> labels_;
  std::vector<std::string> sequences_;
  std::shared_ptr<const MoveTester> move_tester_;
  bool optimize_models_;

  pll::Partition partition_;
  Position default_position_;
  std::queue<Position> starting_positions_;
  std::mutex mutex_;

  // we use a deque here because using a vector would require that
  // wanderers were copy-constructible, which they're not because
  // their partitions can't be copied
  std::deque<Wanderer> wanderers_;
  std::deque<std::future<void>> futures_;
  std::atomic<size_t> idle_wanderer_count_;
  std::vector<bool> wanderer_ready_;

 public:
  Guru(const Options& options,
       pll_utree_t* starting_tree,
       const pll::Model& model,
       const std::vector<std::string>& labels,
       const std::vector<std::string>& sequences);

  Guru(const Options& options,
       const std::vector<pll_utree_t*>& starting_trees,
       const pll::Model& model,
       const std::vector<std::string>& labels,
       const std::vector<std::string>& sequences);

  ~Guru() override;

  bool RequestMove(const Position& position, pll_unode_t* node, MoveType type) override;

  void Start();
  void Wait();

 private:
  void AddStartingPosition(const Position& position);

  std::vector<pll_utree_t*> SortStartingTrees(
      const std::vector<pll_utree_t*>& starting_trees);

  void UpdateWandererStatus(size_t index);
};

} // namespace pt

#endif /* PT_GURU_HPP_ */
