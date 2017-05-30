#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <algorithm>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

#include <model_parameters.hpp>
#include <pll_partition.hpp>
#include <pll_util.hpp>

#include "authority.hpp"
#include "guru.hpp"
#include "wanderer.hpp"

//
// free functions
//

std::vector<std::string> ssplit(const std::string &s, char delim) {
  std::stringstream ss(s);
  std::string item;
  std::vector<std::string> tokens;
  while (std::getline(ss, item, delim)) {
    tokens.push_back(item);
  }
  return tokens;
}

std::map<std::string, double> ReadRaxmlTestData(const std::string& filename)
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

void RunRaxmlTest(const std::string& raxml_data_filename, pt::TreeTable& good_trees)
{
  std::map<std::string, double> raxml_lnls = ReadRaxmlTestData(raxml_data_filename);

  REQUIRE(raxml_lnls.size() == good_trees.size());

  std::vector<std::pair<std::string, double>> raxml_pairs;
  std::vector<std::pair<std::string, double>> pt_pairs;

  for (auto iter = raxml_lnls.begin(); iter != raxml_lnls.end(); ++iter) {
    std::string tree_str = iter->first;
    double raxml_lnl = iter->second;
    double pt_lnl = good_trees.find(tree_str);

    // Test that pt's log-likelihoods are close to RAxML's.
    CHECK(pt_lnl == Approx(raxml_lnl).epsilon(5e-4));

    // Store the (tree_str, lnl) pairs for comparing relative order below.
    raxml_pairs.push_back(std::make_pair(tree_str, raxml_lnl));
    pt_pairs.push_back(std::make_pair(tree_str, pt_lnl));
  }

  // Sort the pt and RAxML (tree_str, lnl) pairs by log-likelihood.
  auto pair_cmp = [](const std::pair<std::string, double>& lhs,
                     const std::pair<std::string, double>& rhs) {
    return lhs.second < rhs.second;
  };

  std::sort(raxml_pairs.begin(), raxml_pairs.end(), pair_cmp);
  std::sort(pt_pairs.begin(), pt_pairs.end(), pair_cmp);

  // Test that the tree ordering is the same for pt and RAxML.
  for (size_t i = 0; i < raxml_pairs.size(); ++i) {
    REQUIRE(raxml_pairs[i].first == pt_pairs[i].first);
  }
}

//
// tests
//

TEST_CASE("partition operations are correct", "[partition]")
{
  std::string newick_path("test-data/tiny/newton.tre");
  std::string fasta_path("test-data/tiny/newton.fasta");
  std::string raxml_path("test-data/tiny/RAxML_info.newton");

  unsigned int tip_node_count;
  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str(),
                                             &tip_node_count);

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, tip_node_count, labels, sequences);

  pt::pll::ModelParameters parameters = pt::pll::ParseRaxmlInfo(raxml_path);

  pt::pll::Partition partition(tree, tip_node_count, parameters, labels, sequences);
  partition.TraversalUpdate(tree, pt::pll::TraversalType::FULL);

  SECTION("log-likelihoods are computed correctly")
  {
    REQUIRE(partition.LogLikelihood(tree) == Approx(-33.387713));
  }

  SECTION("branch lengths are optimized correctly")
  {
    double original_length = tree->length;
    double optimized_length = partition.OptimizeBranch(tree);

    // verify that the length changed
    REQUIRE(original_length != Approx(optimized_length));

    // verify that the length is correct
    REQUIRE(optimized_length == Approx(2.607098));

    // verify that the tree was modified
    REQUIRE(tree->length == optimized_length);
    REQUIRE(tree->back->length == optimized_length);
  }

  SECTION("partitions initialized with a different node order are equivalent")
  {
    // clone the tree in such a way that the tips returned by
    // pll_utree_query_tipnodes() are in a different order
    pll_utree_t* other_tree = pll_utree_clone(tree->back->next->next);
    pll_utree_every(other_tree, pt::pll::cb_copy_clv_traversal);

    pt::pll::Partition other_partition(other_tree, tip_node_count, parameters,
                                       labels, sequences);
    other_partition.TraversalUpdate(other_tree, pt::pll::TraversalType::FULL);

    REQUIRE(other_partition.LogLikelihood(other_tree) ==
            Approx(partition.LogLikelihood(tree)));

    // We did a TraversalUpdate() on this tree, so free its node data.
    pll_utree_every(other_tree, pt::pll::cb_erase_data);
    pll_utree_destroy(other_tree);
  }

  // We did a TraversalUpdate() on this tree, so free its node data.
  pll_utree_every(tree, pt::pll::cb_erase_data);
  pll_utree_destroy(tree);
}

TEST_CASE("wanderer operations are correct", "[wanderer]") {
  std::string newick_path("test-data/five/RAxML_bestTree.five");
  std::string fasta_path("test-data/five/five.fasta");
  std::string raxml_path("test-data/five/RAxML_info.five");

  unsigned int tip_node_count;
  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str(),
                                             &tip_node_count);

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, tip_node_count, labels, sequences);

  pt::pll::ModelParameters parameters = pt::pll::ParseRaxmlInfo(raxml_path);

  pt::pll::Partition partition(tree, tip_node_count, parameters, labels, sequences);
  partition.TraversalUpdate(tree, pt::pll::TraversalType::FULL);

  double ml_lnl = partition.LogLikelihood(tree);

  SECTION("the initial tree's log-likelihood is computed correctly") {
    REQUIRE(ml_lnl == Approx(-3737.47));
  }

  SECTION("wanderers agree with other methods") {
    // from test_pt.cpp: good trees are those with a log-likelihood of
    // at least -3820 (ML is -3737.47).
    double lnl_threshold = -3820.0;
    double lnl_offset = lnl_threshold - ml_lnl;  // -82.53

    pt::Authority authority(ml_lnl, lnl_offset);

    SECTION("using partition move constructor") {
      pt::Wanderer wanderer(authority, std::move(partition), tree);

      wanderer.Start();

      pt::TreeTable& good_trees = authority.GetGoodTreeTable();
      pt::TreeTable& visited_trees = authority.GetVisitedTreeTable();

      SECTION("wanderers agree with old pt") {
        CHECK(good_trees.size() == 13);
        CHECK(visited_trees.size() == 15);

        CHECK(good_trees.contains("(Ref.A1.AU.03.PS1044_Day0.DQ676872,"
                                  "((Ref.A1.RW.92.92RW008.AB253421,"
                                  "Ref.A1.UG.92.92UG037.AB253429),"
                                  "Ref.A2.CM.01.01CM_1445MV.GU201516),"
                                  "Ref.A2.CD.97.97CDKTB48.AF286238);"));
      }

      SECTION("wanderers agree with RAxML") {
        RunRaxmlTest("test-data/five/good_trees.five.raxml", good_trees);
      }

      SECTION("good tree tables are filtered correctly") {
        double new_lnl_threshold = -3800.0;

        // There are 3 good trees with a log-likelihood greater than -3800.
        authority.FilterGoodTreeTable(new_lnl_threshold);

        // We already have a reference to the live table.
        CHECK(good_trees.size() == 3);
      }
    }

    SECTION("using in-place partition constructor") {
      pt::Wanderer wanderer(authority, tree, tip_node_count, parameters,
                            labels, sequences);

      wanderer.Start();

      pt::TreeTable& good_trees = authority.GetGoodTreeTable();
      pt::TreeTable& visited_trees = authority.GetVisitedTreeTable();

      SECTION("wanderers agree with old pt") {
        CHECK(good_trees.size() == 13);
        CHECK(visited_trees.size() == 15);

        CHECK(good_trees.contains("(Ref.A1.AU.03.PS1044_Day0.DQ676872,"
                                  "((Ref.A1.RW.92.92RW008.AB253421,"
                                  "Ref.A1.UG.92.92UG037.AB253429),"
                                  "Ref.A2.CM.01.01CM_1445MV.GU201516),"
                                  "Ref.A2.CD.97.97CDKTB48.AF286238);"));
      }

      SECTION("wanderers agree with RAxML") {
        RunRaxmlTest("test-data/five/good_trees.five.raxml", good_trees);
      }

      SECTION("good tree tables are filtered correctly") {
        double new_lnl_threshold = -3800.0;

        // There are 3 good trees with a log-likelihood greater than -3800.
        authority.FilterGoodTreeTable(new_lnl_threshold);

        // We already have a reference to the live table.
        CHECK(good_trees.size() == 3);
      }
    }
  }

  // We did a TraversalUpdate() on this tree, so free its node data.
  pll_utree_every(tree, pt::pll::cb_erase_data);
  pll_utree_destroy(tree);
}

TEST_CASE("simple guru operations are correct", "[guru_simple]") {
  std::string newick_path("test-data/five/RAxML_bestTree.five");
  std::string fasta_path("test-data/five/five.fasta");
  std::string raxml_path("test-data/five/RAxML_info.five");

  unsigned int tip_node_count;
  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str(),
                                             &tip_node_count);

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, tip_node_count, labels, sequences);

  pt::pll::ModelParameters parameters = pt::pll::ParseRaxmlInfo(raxml_path);

  // from test_pt.cpp: good trees are those with a log-likelihood of
  // at least -3820 (ML is -3737.47).
  double lnl_offset = -82.53;

  SECTION("single-threaded operation is correct") {
    size_t thread_count = 1;

    pt::Guru guru(lnl_offset, thread_count, tree, tip_node_count, parameters,
                  labels, sequences);

    guru.Start();
    guru.Wait();

    pt::TreeTable& good_trees = guru.GetGoodTreeTable();
    pt::TreeTable& visited_trees = guru.GetVisitedTreeTable();

    CHECK(good_trees.size() == 13);
    CHECK(visited_trees.size() == 15);

    CHECK(good_trees.contains("(Ref.A1.AU.03.PS1044_Day0.DQ676872,"
                              "((Ref.A1.RW.92.92RW008.AB253421,"
                              "Ref.A1.UG.92.92UG037.AB253429),"
                              "Ref.A2.CM.01.01CM_1445MV.GU201516),"
                              "Ref.A2.CD.97.97CDKTB48.AF286238);"));
  }

  SECTION("duplicate starting trees don't affect results") {
    size_t thread_count = 2;

    pt::Guru guru(lnl_offset, thread_count, tree, tip_node_count, parameters,
                  labels, sequences);

    // Add the starting tree again. Since we've requested multiple
    // threads, the guru will create wanderers which will go idle
    // immediately upon seeing that their starting tree has already
    // been visited. This shouldn't affect the results.
    guru.AddStartingTree(tree);

    guru.Start();
    guru.Wait();

    pt::TreeTable& good_trees = guru.GetGoodTreeTable();
    pt::TreeTable& visited_trees = guru.GetVisitedTreeTable();

    CHECK(good_trees.size() == 13);
    CHECK(visited_trees.size() == 15);

    CHECK(good_trees.contains("(Ref.A1.AU.03.PS1044_Day0.DQ676872,"
                              "((Ref.A1.RW.92.92RW008.AB253421,"
                              "Ref.A1.UG.92.92UG037.AB253429),"
                              "Ref.A2.CM.01.01CM_1445MV.GU201516),"
                              "Ref.A2.CD.97.97CDKTB48.AF286238);"));
  }

  // Unlike in earlier tests, TraversalUpdate() is never called on
  // this tree, so we don't have to free any node data.
  pll_utree_destroy(tree);
}

TEST_CASE("guru operations on DS1 are correct", "[guru_DS1]") {
  std::string newick_path("test-data/hohna_datasets_fasta/RAxML_bestTree.DS1");
  std::string fasta_path("test-data/hohna_datasets_fasta/DS1.fasta");
  std::string raxml_path("test-data/hohna_datasets_fasta/RAxML_info.DS1");

  unsigned int tip_node_count;
  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str(),
                                             &tip_node_count);

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, tip_node_count, labels, sequences);

  pt::pll::ModelParameters parameters = pt::pll::ParseRaxmlInfo(raxml_path);

  double lnl_offset = -2.0;

  SECTION("single-threaded operation is correct") {
    size_t thread_count = 1;

    pt::Guru guru(lnl_offset, thread_count, tree, tip_node_count, parameters,
                  labels, sequences);

    guru.Start();
    guru.Wait();

    pt::TreeTable& good_trees = guru.GetGoodTreeTable();
    pt::TreeTable& visited_trees = guru.GetVisitedTreeTable();

    CHECK(good_trees.size() == 15);
    CHECK(visited_trees.size() == 659);
  }

  SECTION("multi-threaded operation is correct") {
    size_t thread_count = 4;

    pt::Guru guru(lnl_offset, thread_count, tree, tip_node_count, parameters,
                  labels, sequences);

    guru.Start();
    guru.Wait();

    pt::TreeTable& good_trees = guru.GetGoodTreeTable();
    pt::TreeTable& visited_trees = guru.GetVisitedTreeTable();

    CHECK(good_trees.size() == 15);
    CHECK(visited_trees.size() == 659);
  }

  pll_utree_destroy(tree);
}
