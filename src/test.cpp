// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "partition.hpp"
#include <thread>
namespace pt {
TEST_CASE("Partition", "[partition]") {
  auto p_newton = std::unique_ptr<pt::Partition>(new pt::Partition(
      "test-data/tiny/newton.tre", "test-data/tiny/newton.fasta",
      "test-data/tiny/RAxML_info.newton"));
  double logl = p_newton->FullTraversalLogLikelihood(p_newton->tree_);
  // Value reported by running newton PLL example.
  REQUIRE(fabs(-33.387713 - logl) < 1e-6);

  // Value reported by running newton PLL example.
  REQUIRE(fabs(p_newton->OptimizeCurrentBranch(p_newton->tree_) - 2.607098) <
          1e-6);

  // Verify ToNewick functionality pre-ordering and post-ordering.
  REQUIRE(p_newton->ToNewick(p_newton->tree_) == "((0,1),2,3);");

  // Order and check status.
  p_newton->tree_ = p_newton->ToOrderedNewick(p_newton->tree_);
  REQUIRE(p_newton->ToNewick(p_newton->tree_) == "(0,1,(2,3));");
  p_newton->FullBranchOpt(p_newton->tree_);
  logl = p_newton->FullTraversalLogLikelihood(p_newton->tree_);
  // Tree topologies and likelihoods for all possible NNI moves for newton tree.
  p_newton->MakeTables(2, logl, p_newton->tree_);

  // Print All Tables
  p_newton->PrintTables(1);

  // Check that good tree is in newton good table.
  REQUIRE(p_newton->good_.contains("(0,(1,3),2);"));

  std::cout << "DONE WITH TEST 1" << std::endl;
}

TEST_CASE("NNI_Recursion", "[nnirecursion]") {
  auto p_DS1 = std::unique_ptr<pt::Partition>(
      new pt::Partition("test-data/hohna_datasets_fasta/RAxML_bestTree.DS1",
                        "test-data/hohna_datasets_fasta/DS1.fasta",
                        "test-data/hohna_datasets_fasta/RAxML_info.DS1"));
  // Optimize initial topology.
  p_DS1->FullBranchOpt(p_DS1->tree_);

  // Set ML parameter.
  double logl = p_DS1->FullTraversalLogLikelihood(p_DS1->tree_);

  // Return all trees with a log likelihood of at least -6490 (ML is -6486.9).
  p_DS1->MakeTables(1.000477886, logl, p_DS1->tree_);

  // Only print the good table.
  p_DS1->PrintTables(0);
}
}
