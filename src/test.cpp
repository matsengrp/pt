// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "ctpl_stl.h"
#include "partition.hpp"
namespace pt {
TEST_CASE("Partition", "[partition]") {
  auto p_newton = std::unique_ptr<pt::Partition>(new pt::Partition(
      "test-data/tiny/newton.tre", "test-data/tiny/newton.fasta",
      "test-data/tiny/RAxML_info.newton"));
  TreeTable good_;
  TreeTable all_;
  ctpl::thread_pool pool_(2);
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
  p_newton->MakeTables(4, logl, p_newton->tree_, good_, all_, pool_);

  pool_.stop(true);
  // Print All Tables
  p_newton->PrintTables(1, good_, all_);

  // Check that good tree is in newton good table.
  REQUIRE(good_.contains("(0,(1,3),2);"));

  std::cout << "DONE WITH TEST 1" << std::endl;
  p_newton.reset();
}

TEST_CASE("MultiThreading", "[multithreading]") {

  auto p_five = std::unique_ptr<pt::Partition>(new pt::Partition(
      "test-data/five/RAxML_bestTree.five", "test-data/five/five.fasta",
      "test-data/five/RAxML_info.five"));
  TreeTable good_;
  TreeTable all_;
  ctpl::thread_pool pool_(10);
  // Optimize initial topology.
  p_five->FullBranchOpt(p_five->tree_);

  // Set ML parameter.
  double logl = p_five->FullTraversalLogLikelihood(p_five->tree_);

  // Good trees are trees with a log likelihood of at least -3820 (ML is
  // -3737.47).
  p_five->MakeTables(1.022081783, logl, p_five->tree_, good_, all_, pool_);

  // Wait until all threads in the pool have executed.
  pool_.stop(true);

  // Print both tables.
  // Note that program returns all 15 topologies for a 5-leaf tree.
  p_five->PrintTables(1, good_, all_);

  // Ensure that 13 trees were sorted into the good table, all 15 topologies
  // were investigated.
  REQUIRE(good_.size() == 13);
  REQUIRE(all_.size() == 15);

  // Ensure a specific good tree is contained in good table.
  REQUIRE(good_.contains("(Ref.A1.AU.03.PS1044_Day0.DQ676872,((Ref.A1.RW.92."
                         "92RW008.AB253421,Ref.A1.UG.92.92UG037.AB253429),Ref."
                         "A2.CM.01.01CM_1445MV.GU201516),Ref.A2.CD.97."
                         "97CDKTB48.AF286238);"));

  p_five.reset();
}
TEST_CASE("RAxML info", "[RAxMLinfo]") {
  auto p_five = std::unique_ptr<pt::Partition>(new pt::Partition(
      "test-data/five/RAxML_bestTree.five", "test-data/five/five.fasta",
      "test-data/five/RAxML_info.five"));
  // Values from RAxML info file.
  double rates[6] = {1.420268, 5.241274, 0.870287,
                     0.348533, 7.363575, 1.000000};
  // Verify rates are parsed correctly.
  for (unsigned int i = 0; i < 6; i++) {
    REQUIRE(rates[i] == p_five->GetPartition()->subst_params[0][i]);
  }
  // Values from RAxML info file.
  double frequencies[4] = {.370, 0.194, 0.246, 0.191};
  // Verify rates are parsed correctly.
  for (unsigned int i = 0; i < 4; i++) {
    REQUIRE(frequencies[i] == p_five->GetPartition()->frequencies[0][i]);
  }
}
TEST_CASE("Copy", "[copy]") {
  auto p_five = std::unique_ptr<pt::Partition>(new pt::Partition(
      "test-data/five/RAxML_bestTree.five", "test-data/five/five.fasta",
      "test-data/five/RAxML_info.five"));
  // Copy using the copy constructor, with above partition and its root
  // tree as args.
  auto p_five1 =
      std::unique_ptr<pt::Partition>(new pt::Partition(*p_five, p_five->tree_));
  // Fully optimize both topologies.
  p_five1->FullBranchOpt(p_five1->tree_);
  p_five->FullBranchOpt(p_five->tree_);
  // Check the copy logl is the same as original logl after both have been
  // optimized.
  REQUIRE(fabs(p_five1->FullTraversalLogLikelihood(p_five1->tree_) -
               p_five->FullTraversalLogLikelihood(p_five->tree_)) < 1e-6);
}
}
