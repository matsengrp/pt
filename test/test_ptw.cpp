#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <algorithm>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include <libpll/pll.h>
#include <model.hpp>
#include <pll_partition.hpp>
#include <pll_util.hpp>

#include "authority.hpp"
#include "compressed_tree.hpp"
#include "guru.hpp"
#include "label_dictionary.hpp"
#include "move_tester/always.hpp"
#include "move_tester/branch_neighborhood_optimizer.hpp"
#include "move_tester/single_branch_optimizer.hpp"
#include "options.hpp"
#include "ordered_tree.hpp"
#include "position.hpp"
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

double FindByString(pt::TreeTable& trees, const std::string& tree_str)
{
  auto table = trees.lock_table();
  for (auto& item : table) {
    if (item.first.Decode() == tree_str) {
      return item.second;
    }
  }

  throw std::invalid_argument("tree string not found in table");
}

double ContainsByString(pt::TreeTable& trees, const std::string& tree_str)
{
  auto table = trees.lock_table();
  for (auto& item : table) {
    if (item.first.Decode() == tree_str) {
      return true;
    }
  }

  return false;
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
    double pt_lnl = FindByString(good_trees, tree_str);

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

TEST_CASE("wanderer operations are correct", "[wanderer]") {
  std::string newick_path("test-data/five/RAxML_bestTree.five");
  std::string fasta_path("test-data/five/five.fasta");
  std::string raxml_path("test-data/five/RAxML_info.five");

  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str());

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, tree->tip_count, labels, sequences);

  pt::CompressedTree::BuildDictionary(labels);

  pt::pll::Model model = pt::pll::ParseRaxmlInfo(raxml_path);

  pt::pll::Partition partition(tree, model, labels, sequences);

  pll_unode_t* node = tree->nodes[tree->tip_count + tree->inner_count - 1];
  partition.TraversalUpdate(node, pt::pll::TraversalType::FULL);

  double ml_lnl = partition.LogLikelihood(node);

  SECTION("the initial tree's log-likelihood is computed correctly") {
    REQUIRE(ml_lnl == Approx(-3737.47));
  }

  SECTION("wanderers agree with other methods") {
    // from test_pt.cpp: good trees are those with a log-likelihood of
    // at least -3820 (ML is -3737.47).
    double lnl_threshold = -3820.0;
    double lnl_offset = lnl_threshold - ml_lnl;  // -82.53

    pt::Options options;
    options.lnl_offset = lnl_offset;
    options.move_tester = std::make_shared<pt::move_tester::Always>();
    options.optimize_models = false;

    pt::Authority authority(options, ml_lnl);

    SECTION("using in-place partition constructor") {
      pt::Wanderer wanderer(authority, pt::Position(tree, model), labels, sequences,
                            options.move_tester, options.optimize_models);

      wanderer.Start();

      pt::TreeTable& good_trees = authority.GetGoodTreeTable();
      pt::TreeTable& visited_trees = authority.GetVisitedTreeTable();
      pt::TreeTable& tested_trees = authority.GetTestedTreeTable();

      SECTION("wanderers agree with old pt") {
        CHECK(good_trees.size() == 13);
        CHECK(visited_trees.size() == 15);
        CHECK(tested_trees.size() == 14);

        CHECK(ContainsByString(good_trees,
                               "(Ref.A1.AU.03.PS1044_Day0.DQ676872,"
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
  pll_utree_destroy(tree, pt::pll::cb_erase_data);
}

TEST_CASE("simple guru operations are correct", "[guru_simple]") {
  std::string newick_path("test-data/five/RAxML_bestTree.five");
  std::string fasta_path("test-data/five/five.fasta");
  std::string raxml_path("test-data/five/RAxML_info.five");

  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str());

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, tree->tip_count, labels, sequences);

  pt::CompressedTree::BuildDictionary(labels);

  pt::pll::Model model = pt::pll::ParseRaxmlInfo(raxml_path);

  // from test_pt.cpp: good trees are those with a log-likelihood of
  // at least -3820 (ML is -3737.47).
  double lnl_offset = -82.53;

  pt::Options options;
  options.lnl_offset = lnl_offset;
  options.move_tester = std::make_shared<pt::move_tester::Always>();

  SECTION("single-threaded operation is correct") {
    options.thread_count = 1;

    pt::Guru guru(options, tree, model, labels, sequences);

    guru.Start();
    guru.Wait();

    pt::TreeTable& good_trees = guru.GetGoodTreeTable();
    pt::TreeTable& visited_trees = guru.GetVisitedTreeTable();
    pt::TreeTable& tested_trees = guru.GetTestedTreeTable();

    CHECK(good_trees.size() == 13);
    CHECK(visited_trees.size() == 15);
    CHECK(tested_trees.size() == 14);

    CHECK(ContainsByString(good_trees,
                           "(Ref.A1.AU.03.PS1044_Day0.DQ676872,"
                           "((Ref.A1.RW.92.92RW008.AB253421,"
                           "Ref.A1.UG.92.92UG037.AB253429),"
                           "Ref.A2.CM.01.01CM_1445MV.GU201516),"
                           "Ref.A2.CD.97.97CDKTB48.AF286238);"));
  }

  SECTION("duplicate starting trees don't affect results") {
    options.thread_count = 2;

    // Populate a vector with duplicate starting trees.
    std::vector<pll_utree_t*> trees(2, tree);

    pt::Guru guru(options, trees, model, labels, sequences);

    // Since we've requested multiple threads, the guru will create
    // wanderers which will go idle immediately upon seeing that their
    // starting tree has already been visited. This shouldn't affect
    // the results.

    guru.Start();
    guru.Wait();

    pt::TreeTable& good_trees = guru.GetGoodTreeTable();
    pt::TreeTable& visited_trees = guru.GetVisitedTreeTable();
    pt::TreeTable& tested_trees = guru.GetTestedTreeTable();

    CHECK(good_trees.size() == 13);
    CHECK(visited_trees.size() == 15);
    CHECK(tested_trees.size() == 14);

    CHECK(ContainsByString(good_trees,
                           "(Ref.A1.AU.03.PS1044_Day0.DQ676872,"
                           "((Ref.A1.RW.92.92RW008.AB253421,"
                           "Ref.A1.UG.92.92UG037.AB253429),"
                           "Ref.A2.CM.01.01CM_1445MV.GU201516),"
                           "Ref.A2.CD.97.97CDKTB48.AF286238);"));
  }

  // Unlike in earlier tests, TraversalUpdate() is never called on
  // this tree, so we don't have to free any node data.
  pll_utree_destroy(tree, nullptr);
}

void RunGuruTest(const pt::Options& options,
                 pll_utree_t* tree,
                 const pt::pll::Model& model,
                 const std::vector<std::string>& labels,
                 const std::vector<std::string>& sequences,
                 size_t good_tree_count,
                 size_t visited_tree_count,
                 size_t tested_tree_count)
{
  pt::Guru guru(options, tree, model, labels, sequences);

  guru.Start();
  guru.Wait();

  pt::TreeTable& good_trees = guru.GetGoodTreeTable();
  pt::TreeTable& visited_trees = guru.GetVisitedTreeTable();
  pt::TreeTable& tested_trees = guru.GetTestedTreeTable();

  CHECK(good_trees.size() == good_tree_count);
  CHECK(visited_trees.size() == visited_tree_count);
  CHECK(tested_trees.size() == tested_tree_count);
}

TEST_CASE("guru operations on DS1 are correct", "[guru_DS1]") {
  std::string newick_path("test-data/hohna_datasets_fasta/RAxML_bestTree.DS1");
  std::string fasta_path("test-data/hohna_datasets_fasta/DS1.fasta");
  std::string raxml_path("test-data/hohna_datasets_fasta/RAxML_info.DS1");

  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str());

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, tree->tip_count, labels, sequences);

  pt::CompressedTree::BuildDictionary(labels);

  pt::pll::Model model = pt::pll::ParseRaxmlInfo(raxml_path);

  pt::Options options;

  SECTION("using move_tester::Always") {
    options.lnl_offset = -2.0;
    options.move_tester = std::make_shared<pt::move_tester::Always>();

    size_t good_tree_count = 15;
    size_t visited_tree_count = 659;
    size_t tested_tree_count = 658;

    SECTION("single-threaded operation is correct") {
      options.thread_count = 1;

      RunGuruTest(options, tree, model, labels, sequences,
                  good_tree_count, visited_tree_count, tested_tree_count);
    }

    SECTION("multi-threaded operation is correct") {
      options.thread_count = 4;

      RunGuruTest(options, tree, model, labels, sequences,
                  good_tree_count, visited_tree_count, tested_tree_count);
    }
  }

  SECTION("using move_tester::SingleBranchOptimizer") {
    options.lnl_offset = -2.0;
    options.move_tester = std::make_shared<pt::move_tester::SingleBranchOptimizer>();

    size_t good_tree_count = 9;
    size_t visited_tree_count = 9;
    size_t tested_tree_count = 404;

    SECTION("single-threaded operation is correct") {
      options.thread_count = 1;

      RunGuruTest(options, tree, model, labels, sequences,
                  good_tree_count, visited_tree_count, tested_tree_count);
    }

    SECTION("multi-threaded operation is correct") {
      options.thread_count = 4;

      RunGuruTest(options, tree, model, labels, sequences,
                  good_tree_count, visited_tree_count, tested_tree_count);
    }
  }

  SECTION("using move_tester::BranchNeighborhoodOptimizer") {
    options.lnl_offset = -2.0;

    SECTION("a non-positive radius is not accepted") {
      CHECK_THROWS_AS(pt::move_tester::BranchNeighborhoodOptimizer tmp(-1),
                      std::invalid_argument);
      CHECK_THROWS_AS(pt::move_tester::BranchNeighborhoodOptimizer tmp(0),
                      std::invalid_argument);
    }

    SECTION("with an optimization radius of 1") {
      int optimization_radius = 1;
      options.move_tester =
          std::make_shared<pt::move_tester::BranchNeighborhoodOptimizer>(optimization_radius);

      size_t good_tree_count = 14;
      size_t visited_tree_count = 14;
      size_t tested_tree_count = 617;

      SECTION("single-threaded operation is correct") {
        options.thread_count = 1;

        RunGuruTest(options, tree, model, labels, sequences,
                    good_tree_count, visited_tree_count, tested_tree_count);
      }

      SECTION("multi-threaded operation is correct") {
        options.thread_count = 4;

        RunGuruTest(options, tree, model, labels, sequences,
                    good_tree_count, visited_tree_count, tested_tree_count);
      }
    }

    SECTION("with an optimization radius of 2") {
      int optimization_radius = 2;
      options.move_tester =
          std::make_shared<pt::move_tester::BranchNeighborhoodOptimizer>(optimization_radius);

      size_t good_tree_count = 15;
      size_t visited_tree_count = 15;
      size_t tested_tree_count = 658;

      SECTION("single-threaded operation is correct") {
        options.thread_count = 1;

        RunGuruTest(options, tree, model, labels, sequences,
                    good_tree_count, visited_tree_count, tested_tree_count);
      }

      SECTION("multi-threaded operation is correct") {
        options.thread_count = 4;

        RunGuruTest(options, tree, model, labels, sequences,
                    good_tree_count, visited_tree_count, tested_tree_count);
      }
    }
  }

  SECTION("when tested tree tracking is disabled") {
    options.lnl_offset = -2.0;
    options.move_tester = std::make_shared<pt::move_tester::Always>();
    options.track_tested_trees = false;

    size_t good_tree_count = 15;
    size_t visited_tree_count = 659;
    size_t tested_tree_count = 0;

    SECTION("single-threaded operation is correct") {
      options.thread_count = 1;

      RunGuruTest(options, tree, model, labels, sequences,
                  good_tree_count, visited_tree_count, tested_tree_count);
    }

    SECTION("multi-threaded operation is correct") {
      options.thread_count = 4;

      RunGuruTest(options, tree, model, labels, sequences,
                  good_tree_count, visited_tree_count, tested_tree_count);
    }
  }

  SECTION("with only one rate category") {
    options.lnl_offset = -2.0;
    options.move_tester = std::make_shared<pt::move_tester::SingleBranchOptimizer>();
    options.thread_count = 1;
    options.rate_categories = 1;

    model = pt::pll::ParseRaxmlInfo(raxml_path, options.rate_categories);

    size_t good_tree_count = 5;
    size_t visited_tree_count = 8;
    size_t tested_tree_count = 230;

    RunGuruTest(options, tree, model, labels, sequences,
                good_tree_count, visited_tree_count, tested_tree_count);
  }

  SECTION("with model optimization enabled") {
    options.lnl_offset = -3.0;
    options.move_tester = std::make_shared<pt::move_tester::SingleBranchOptimizer>();
    options.thread_count = 1;
    options.optimize_models = true;

    size_t good_tree_count = 45;
    size_t visited_tree_count = 45;
    size_t tested_tree_count = 1873;

    RunGuruTest(options, tree, model, labels, sequences,
                good_tree_count, visited_tree_count, tested_tree_count);
  }

  pll_utree_destroy(tree, nullptr);
}

TEST_CASE("guru operations on DS1 with two starting trees are correct", "[guru_DS1_two_trees]") {
  std::string newick_path("test-data/hohna_datasets_fasta/two-peaks/two_peaks.nw");

  std::string fasta_path("test-data/hohna_datasets_fasta/DS1.fasta");
  std::string raxml_path("test-data/hohna_datasets_fasta/RAxML_info.DS1");

  std::vector<pll_utree_t*> trees = pt::pll::ParseMultiNewick(newick_path);

  REQUIRE(trees.size() == 2);

  pll_utree_t* first_tree = trees[0];
  pll_utree_t* second_tree = trees[1];

  REQUIRE(first_tree);
  REQUIRE(second_tree);

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(fasta_path, first_tree->tip_count, labels, sequences);

  pt::CompressedTree::BuildDictionary(labels);

  pt::pll::Model model = pt::pll::ParseRaxmlInfo(raxml_path);

  pt::Options options;
  options.thread_count = 1;
  options.move_tester = std::make_shared<pt::move_tester::Always>();

  SECTION("on the first peak") {
    options.lnl_offset = -2.0;

    pt::Guru guru(options, first_tree, model, labels, sequences);

    guru.Start();
    guru.Wait();

    pt::TreeTable& good_trees = guru.GetGoodTreeTable();
    pt::TreeTable& visited_trees = guru.GetVisitedTreeTable();

    CHECK(good_trees.size() == 8);
    CHECK(visited_trees.size() == 353);
  }

  SECTION("on the second peak") {
    options.lnl_offset = -1.05;

    pt::Guru guru(options, second_tree, model, labels, sequences);

    guru.Start();
    guru.Wait();

    pt::TreeTable& good_trees = guru.GetGoodTreeTable();
    pt::TreeTable& visited_trees = guru.GetVisitedTreeTable();

    CHECK(good_trees.size() == 4);
    CHECK(visited_trees.size() == 185);
  }

  SECTION("on both peaks, given the higher peak first") {
    options.lnl_offset = -2.0;

    pt::Guru guru(options, trees, model, labels, sequences);

    guru.Start();
    guru.Wait();

    pt::TreeTable& good_trees = guru.GetGoodTreeTable();
    pt::TreeTable& visited_trees = guru.GetVisitedTreeTable();

    CHECK(good_trees.size() == 12);
    CHECK(visited_trees.size() == 538);
  }

  SECTION("on both peaks, given the lower peak first") {
    options.lnl_offset = -2.0;

    std::vector<pll_utree_t*> reordered_trees{trees[1], trees[0]};

    pt::Guru guru(options, reordered_trees, model, labels, sequences);

    guru.Start();
    guru.Wait();

    pt::TreeTable& good_trees = guru.GetGoodTreeTable();
    pt::TreeTable& visited_trees = guru.GetVisitedTreeTable();

    CHECK(good_trees.size() == 12);
    CHECK(visited_trees.size() == 538);
  }

  for (auto tree : trees) {
    pll_utree_destroy(tree, nullptr);
  }
}

// TODO: this is duplicated in authority.cpp
std::string OrderedNewickString(pll_utree_t* tree)
{
  // clone the tree, because ToOrderedNewick() reorders the tree in
  // place and we don't want to modify the tree we were given. note
  // that we aren't modifying any of the node user data (via the
  // node->data pointers) so we don't have to copy or free those.
  pll_utree_t* clone = pll_utree_clone(tree);

  // ToOrderedNewick() only reorders the tree, despite its name. its
  // return value is the rerooted and reordered tree, which we assign
  // to our pointer before continuing.
  pll_unode_t* root = pt::ToOrderedNewick(pt::GetVirtualRoot(clone));
  std::string newick_str = pt::ToNewick(root);

  pll_utree_destroy(clone, nullptr);
  return newick_str;
}

TEST_CASE("label dictionaries work correctly", "[label_dictionary]") {
  SECTION("in normal cases") {
    std::vector<std::string> labels{"a", "b", "c"};
    pt::LabelDictionary dict(labels);

    CHECK(dict.GetIndex("a") == 0);
    CHECK(dict.GetIndex("b") == 1);
    CHECK(dict.GetIndex("c") == 2);

    CHECK(dict.GetLabel(0) == "a");
    CHECK(dict.GetLabel(1) == "b");
    CHECK(dict.GetLabel(2) == "c");

    CHECK_THROWS_AS(dict.GetIndex("bad"), std::invalid_argument);
    CHECK_THROWS_AS(dict.GetLabel(3), std::invalid_argument);
  }

  SECTION("given duplicate labels") {
    std::vector<std::string> labels{"a", "b", "a"};
    CHECK_THROWS_AS(pt::LabelDictionary dict(labels), std::invalid_argument);
  }

  SECTION("given no labels") {
    std::vector<std::string> labels;
    CHECK_THROWS_AS(pt::LabelDictionary dict(labels), std::invalid_argument);
  }

  SECTION("given too many labels") {
    using index_type = pt::LabelDictionary::index_type;
    size_t max_size = static_cast<size_t>(std::numeric_limits<index_type>::max()) + 1;
    std::vector<std::string> labels(max_size);

    // generate as many unique labels as we need
    size_t i = 0;
    std::generate(labels.begin(), labels.end(), [&i]() {
            std::stringstream ss;
            ss << i++;
            return ss.str();
        });

    CHECK_NOTHROW(pt::LabelDictionary dict(labels));

    labels.push_back("one too many");
    CHECK_THROWS_AS(pt::LabelDictionary dict(labels), std::length_error);
  }
}

TEST_CASE("tree compression works correctly", "[compressed_tree]") {
  std::string newick_path("test-data/five/RAxML_bestTree.five");
  pll_utree_t* tree = pll_utree_parse_newick(newick_path.c_str());

  // build the dictionary
  std::vector<std::string> labels;
  for (size_t i = 0; i < tree->tip_count; ++i) {
    pll_unode_t* node = tree->nodes[i];

    REQUIRE(node->label);
    labels.emplace_back(node->label);
  }

  pt::CompressedTree::BuildDictionary(labels);

  SECTION("given a null tree pointer") {
    CHECK_THROWS_AS(pt::CompressedTree tmp(nullptr),
                    std::invalid_argument);
  }

  SECTION("encoding and decoding works correctly") {
    pt::CompressedTree ct(tree);

    //std::cerr << ct.ToDebugString();

    std::string expected_newick_str = OrderedNewickString(tree);
    std::string actual_newick_str = ct.Decode();

    CHECK(actual_newick_str == expected_newick_str);
  }

  SECTION("different trees are encoded differently") {
    pll_unode_t* root = pt::GetVirtualRoot(tree);

    // get to a non-pendant branch
    root = root->next->back;

    // verify that root and root->back are internal nodes
    REQUIRE(root->next);
    REQUIRE(root->back->next);

    //
    // encode the tree
    //

    pt::CompressedTree ct1(tree);

#if 0
    std::cerr << "original tree ("
              << ct1.Hash() << ")\n\n"
              << ct1.ToDebugString() << "\n";
    pll_utree_show_ascii(root, PLL_UTREE_SHOW_LABEL);
    std::cerr << "\n\n";
#endif

    //
    // apply an NNI move and verify that the encodings differ
    //

    pll_utree_nni(root, pt::MoveType::LEFT, nullptr);

    pt::CompressedTree ct2(tree);

#if 0
    std::cerr << "modified tree ("
              << ct2.Hash() << ")\n\n"
              << ct2.ToDebugString() << "\n";
    pll_utree_show_ascii(root, PLL_UTREE_SHOW_LABEL);
    std::cerr << "\n\n";
#endif

    CHECK(ct1 != ct2);
    CHECK(ct1.Hash() != ct2.Hash());

    //
    // undo the NNI move and verify that we get the original encoding
    //

    pll_utree_nni(root, pt::MoveType::LEFT, nullptr);

    pt::CompressedTree ct3(tree);

#if 0
    std::cerr << "undo tree ("
              << ct3.Hash() << ")\n\n"
              << ct3.ToDebugString() << "\n";
    pll_utree_show_ascii(root, PLL_UTREE_SHOW_LABEL);
    std::cerr << "\n\n";
#endif

    CHECK(ct3 == ct1);
    CHECK(ct3.Hash() == ct1.Hash());

    CHECK(ct3 != ct2);
    CHECK(ct3.Hash() != ct2.Hash());
  }

  pll_utree_destroy(tree, nullptr);
}
