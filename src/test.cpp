// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "partition.hpp"

namespace pt {

TEST_CASE("Partition", "[partition]") {
  auto p_DS1 = std::unique_ptr<Partition>(
      new Partition("test-data/hohna_datasets_fasta/RAxML_bestTree.DS1",
                    "test-data/hohna_datasets_fasta/DS1.fasta",
                    "test-data/hohna_datasets_fasta/RAxML_info.DS1"));

  auto p_newton = std::unique_ptr<Partition>(
      new Partition("test-data/tiny/newton.tre", "test-data/tiny/newton.fasta",
                    "test-data/tiny/RAxML_info.newton"));
  // Value reported by running newton PLL example.
  // REQUIRE(-33.387713 - p_newton->FullTraversalLogLikelihood() < 1e-6);

  // Tree topologies and likelihoods pre and post-NNI.
  p_newton->MakeTables();
}
}
