// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "partition.hpp"
#include <thread>
namespace pt {
TEST_CASE("Partition", "[partition]") {
  auto p_DS1 = std::unique_ptr<pt::Partition>(
      new pt::Partition("test-data/hohna_datasets_fasta/RAxML_bestTree.DS1",
                        "test-data/hohna_datasets_fasta/DS1.fasta",
                        "test-data/hohna_datasets_fasta/RAxML_info.DS1"));

  auto p_newton = std::unique_ptr<pt::Partition>(new pt::Partition(
      "test-data/tiny/newton.tre", "test-data/tiny/newton.fasta",
      "test-data/tiny/RAxML_info.newton"));

  // Value reported by running newton PLL example.
  REQUIRE(fabs(-33.387713 -
               p_newton->FullTraversalLogLikelihood(p_newton->tree_)) < 1e-6);

  // Value reported by running newton PLL example.
  REQUIRE(fabs(p_newton->OptimizeCurrentBranch(p_newton->tree_) - 2.607098) <
          1e-6);

  // Verify ToNewick functionality pre-ordering and post-ordering.
  REQUIRE(p_newton->ToNewick(p_newton->tree_) == "((0,1),2,3);");

  // Order and check status.
  p_newton->tree_ = p_newton->ToOrderedNewick(p_newton->tree_);
  REQUIRE(p_newton->ToNewick(p_newton->tree_) == "(0,1,(2,3));");

  // Tree topologies and likelihoods for all possible NNI moves for newton tree.
  p_newton->MakeTables();
  // Check that good tree is in newton good table.
  // NOTE: Changing the cutoff_ global parameter affects this test.
  REQUIRE(p_newton->good_.contains("(0,(1,3),2);"));
}
}
