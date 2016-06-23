// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "partition.hpp"


namespace pt {

TEST_CASE("Partition", "[partition]") {
  auto p_five = std::unique_ptr<Partition>(new Partition(
      "test-data/five/RAxML_bestTree.five", "test-data/five/five.fasta","test-data/five/RAxML_info.five"));

  REQUIRE(p_five->branch_count() == 7);
  auto p_newton = std::unique_ptr<Partition>(new Partition(
      "test-data/tiny/newton.tre", "test-data/tiny/newton.fasta","test-data/tiny/RAxML_info.newton"));
  // Value reported by running newton PLL example.
  REQUIRE(-33.387713 - p_newton->FullTraversalLogLikelihood() < 1e-6);
}

TEST_CASE("PreOrder", "[preorder]") {
   auto p_newton = std::unique_ptr<Partition>(new Partition(
      "test-data/tiny/newton.tre", "test-data/tiny/newton.fasta","test-data/tiny/RAxML_info.newton"));

      p_newton->BranchLengthsTest1();
}
TEST_CASE("PostOrder", "[postorder]") {
   auto p_newton = std::unique_ptr<Partition>(new Partition(
      "test-data/tiny/newton.tre", "test-data/tiny/newton.fasta","test-data/tiny/RAxML_info.newton"));
      p_newton->ToOrderedNewick();
      p_newton->BranchLengthsTest2();
}

TEST_CASE("TreeNoodle", "[tree]") {
  unsigned int tip_count;
  pll_utree_t* tree =
      pll_utree_parse_newick("test-data/tiny/three.tre", &tip_count);

  ///pll_utree_show_ascii(tree, 255);

  pll_utree_destroy(tree);
}
}

