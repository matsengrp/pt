// This tells Catch to provide a main() - only do this in one cpp file.
#define CATCH_CONFIG_MAIN

#include "catch.hpp"
#include "ctpl_stl.h"
#include "ordered-tree.hpp"
#include "partition.hpp"
namespace pt {
TEST_CASE("Partition", "[partition]") {
  auto p_newton = pt::Partition::Create("test-data/tiny/newton.tre",
                                        "test-data/tiny/newton.fasta",
                                        "test-data/tiny/RAxML_info.newton");
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
  REQUIRE(ToNewick(p_newton->tree_) == "((0,1),2,3);");

  // Order and check status.
  p_newton->tree_ = ToOrderedNewick(p_newton->tree_);
  REQUIRE(ToNewick(p_newton->tree_) == "(0,1,(2,3));");
  p_newton->FullBranchOpt(p_newton->tree_);
  logl = p_newton->FullTraversalLogLikelihood(p_newton->tree_);
  // Tree topologies and likelihoods for all possible NNI moves for newton tree.
  p_newton->MakeTables(4, logl, p_newton->tree_, good_, all_, pool_);

  pool_.stop(true);
  // Print All Tables
  p_newton->PrintTables(1, good_, all_);

  // Check that good tree is in newton good table.
  REQUIRE(good_.contains("(0,(1,3),2);"));
  p_newton.reset();
}

std::map<std::string, double> ReadRaxmlTest(const std::string& filename)
{
  std::map<std::string, double> tree_lnls;

  std::ifstream file(filename);
  std::string line;

  while (std::getline(file, line)) {
    std::vector<std::string> tokens = ssplit(line, '\t');
    tree_lnls[tokens[0]] = std::stod(tokens[1]);
  }

  return tree_lnls;
}

TEST_CASE("MultiThreading", "[multithreading]") {
  auto p_five = pt::Partition::Create("test-data/five/RAxML_bestTree.five",
                                      "test-data/five/five.fasta",
                                      "test-data/five/RAxML_info.five");
  TreeTable good_;
  TreeTable all_;
  ctpl::thread_pool pool_(10);
  // Optimize initial topology.
  p_five->FullBranchOpt(p_five->tree_);

  // Set ML parameter.
  double logl = p_five->FullTraversalLogLikelihood(p_five->tree_);

  REQUIRE(logl == Approx(-3737.47));

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

  // Read in known good trees and their associated RAxML log-likelihoods.
  // See test/test-data/five/README.md for details on this file.
  std::map<std::string, double> raxml_lnls =
      ReadRaxmlTest("test-data/five/good_trees.five.raxml");

  REQUIRE(raxml_lnls.size() == good_.size());

  for (auto iter = raxml_lnls.begin(); iter != raxml_lnls.end(); ++iter) {
    std::string tree = iter->first;
    double raxml_lnl = iter->second;
    double pt_lnl = good_.find(tree);

    // Test that pt's log-likelihoods are close to RAxML's.
    REQUIRE(pt_lnl == Approx(raxml_lnl).epsilon(5e-4));
  }

  p_five.reset();
}
TEST_CASE("RAxML info", "[RAxMLinfo]") {
  auto p_five = pt::Partition::Create("test-data/five/RAxML_bestTree.five",
                                      "test-data/five/five.fasta",
                                      "test-data/five/RAxML_info.five");
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
  p_five.reset();
}
TEST_CASE("Copy", "[copy]") {
  auto p_five = pt::Partition::Create("test-data/five/RAxML_bestTree.five",
                                      "test-data/five/five.fasta",
                                      "test-data/five/RAxML_info.five");
  // Copy using the copy constructor, with above partition and its root
  // tree as args.
  auto p_five1 = pt::Partition::Create(*p_five, p_five->tree_);
  // Fully optimize both topologies.
  p_five1->FullBranchOpt(p_five1->tree_);
  p_five->FullBranchOpt(p_five->tree_);
  // Check the copy logl is the same as original logl after both have been
  // optimized.
  REQUIRE(fabs(p_five1->FullTraversalLogLikelihood(p_five1->tree_) -
               p_five->FullTraversalLogLikelihood(p_five->tree_)) < 1e-6);
  p_five.reset();
}

TEST_CASE("BigExample", "[BigExample]") {
  auto p_DS1 = pt::Partition::Create("test-data/hohna_datasets_fasta/RAxML_bestTree.DS1",
                                     "test-data/hohna_datasets_fasta/DS1.fasta",
                                     "test-data/hohna_datasets_fasta/RAxML_info.DS1");
  auto p_DS1_2 = pt::Partition::Create("test-data/hohna_datasets_fasta/RAxML_secondpeak.DS1",
                                       "test-data/hohna_datasets_fasta/DS1.fasta",
                                       "test-data/hohna_datasets_fasta/RAxML_secondinfo.DS1");
  TreeTable good_;
  TreeTable all_;
  ctpl::thread_pool pool_(6);
  // Optimize initial topology.
  p_DS1->FullBranchOpt(p_DS1->tree_);
  // Set ML parameter. (ML is 6482.13).
  p_DS1->tree_ = ToOrderedNewick(p_DS1->tree_);
  double logl = p_DS1->FullTraversalLogLikelihood(p_DS1->tree_);
  // Add to tables.
  good_.insert(ToNewick(p_DS1->tree_), logl);
  all_.insert(ToNewick(p_DS1->tree_), 0);
  // Explore peak.
  p_DS1->MakeTables(1.000001, logl, p_DS1->tree_, good_, all_, pool_);

  // Repeat for second peak.
  p_DS1_2->FullBranchOpt(p_DS1_2->tree_);
  p_DS1_2->tree_ = ToOrderedNewick(p_DS1_2->tree_);
  good_.insert(ToNewick(p_DS1_2->tree_),
               p_DS1_2->FullTraversalLogLikelihood(p_DS1_2->tree_));
  all_.insert(ToNewick(p_DS1_2->tree_), 0);
  p_DS1_2->MakeTables(1.000001, logl, p_DS1_2->tree_, good_, all_, pool_);
  // Wait until all threads in the pool have executed.
  pool_.stop(true);

  // Print good table.
  p_DS1->PrintTables(0, good_, all_);
  p_DS1.reset();
  p_DS1_2.reset();
}

TEST_CASE("partial likelihoods are evaluated correctly", "[partial]") {
  auto p_five = pt::Partition::Create("test-data/five/RAxML_bestTree.five",
                                      "test-data/five/five.fasta",
                                      "test-data/five/RAxML_info.five");

  // The original length of the edge pointed to by p_five->tree_ is
  // approximately 0.0302113. Set it to a new value so we have
  // something to optimize.
  p_five->tree_->length = 0.1;
  p_five->tree_->back->length = 0.1;

  // Perform the first traversal of the tree. Even though we pass
  // is_full = false here, we'll get a full traversal since none of
  // the nodes have node_info data associated with them yet.
  p_five->TraversalUpdate(p_five->tree_, false);

  // Optimize the current branch.
  double new_length = p_five->OptimizeCurrentBranch(p_five->tree_);

  REQUIRE(new_length == Approx(0.0302113).epsilon(1e-3));

  // Compute the new log-likelihood of the tree as-is.
  double lnl_as_is = p_five->LogLikelihood(p_five->tree_);

  // Recompute the log-likelihood of the tree with a partial traversal.
  p_five->TraversalUpdate(p_five->tree_, false);
  double lnl_partial = p_five->LogLikelihood(p_five->tree_);

  // Recompute the log-likelihood of the tree with a full traversal.
  p_five->TraversalUpdate(p_five->tree_, true);
  double lnl_full = p_five->LogLikelihood(p_five->tree_);

  REQUIRE(lnl_partial == Approx(lnl_full));

  // This isn't true -- a traversal after optimization is required.
  //REQUIRE(lnl_as_is == Approx(lnl_full));
}

} // namespace pt
