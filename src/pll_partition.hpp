#ifndef PT_PLL_PARTITION_HPP_
#define PT_PLL_PARTITION_HPP_

#include <string>
#include <memory>
#include <vector>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

#include "model_parameters.hpp"

namespace pt { namespace pll {

enum class TraversalType { FULL, PARTIAL };


class Partition
{
 private:
  using PartitionPtr = std::unique_ptr<pll_partition_t,
                                       decltype(&pll_partition_destroy)>;

  const unsigned int tip_node_count_;
  PartitionPtr partition_;
  std::vector<unsigned int> params_indices_;

  // scratch buffers for TraversalUpdate()
  std::vector<pll_utree_t*> travbuffer_;
  std::vector<double> branch_lengths_;
  std::vector<unsigned int> matrix_indices_;
  std::vector<pll_operation_t> operations_;

  // scratch buffers for OptimizeBranch()
  double* sumtable_;

 public:
  // based on old Partition constructor in partition.cpp
  Partition(pll_utree_t* tree, unsigned int tip_node_count,
            const ModelParameters& parameters,
            const std::vector<std::string>& labels,
            const std::vector<std::string>& sequences);

  // move constructor
  Partition(Partition&& rhs);

  // explicitly mark other operations as deleted
  Partition(const Partition& rhs) = delete;
  Partition& operator=(const Partition& rhs) = delete;
  Partition& operator=(Partition&& rhs) = delete;

  ~Partition();

  double LogLikelihood(pll_utree_t* tree);

  unsigned int TraversalUpdate(pll_utree_t* root, TraversalType type);
  void UpdateBranchLength(pll_utree_t* node, double length);

  double OptimizeBranch(pll_utree_t* node);
  void OptimizeAllBranchesOnce(pll_utree_t* tree);
  void OptimizeAllBranches(pll_utree_t* tree);

  unsigned int tip_node_count() const;
  unsigned int inner_node_count() const;
  unsigned int node_count() const;
  unsigned int branch_count() const;

 private:
  // based on SetModelParameters() in pll-utils.cpp, but the
  // parameters struct will be populated elsewhere -- no need to parse
  // the RAxML info file here
  void SetModelParameters(const ModelParameters& parameters);

  // based on EquipPartitionWithData() in pll-utils.cpp
  void SetTipStates(pll_utree_t* tree, const std::vector<std::string>& labels,
                    const std::vector<std::string>& sequences);

  void AllocateScratchBuffers();
  void FreeScratchBuffers();
};

//
// Partition inlines
//

inline unsigned int Partition::tip_node_count() const
{
  return tip_node_count_;
}

inline unsigned int Partition::inner_node_count() const
{
  return tip_node_count_ - 2;
}

inline unsigned int Partition::node_count() const
{
  return tip_node_count() + inner_node_count();
}

inline unsigned int Partition::branch_count() const
{
  return node_count() - 1;
}

} } // namespace pt::pll

#endif // PT_PLL_PARTITION_HPP_
