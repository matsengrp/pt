#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <string>

// pll.h is missing a header guard
#ifndef LIBPLL_PLL_H_
#define LIBPLL_PLL_H_
#include <libpll/pll.h>
#endif

#include "model_parameters.hpp"
#include "pll-utils.hpp"
#include "pll_partition.hpp"
#include "wanderer.hpp"

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
  pt::ParseFasta(fasta_path, tip_node_count, labels, sequences);

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

  // FIXME: this is just temporary, to see if a Wanderer can be
  //        created and destroyed
  SECTION("wanderers are created correctly")
  {
    pt::Authority authority(-33.387713, 1.1);

    // FIXME: the wanderer takes ownership of tree here. this is
    //        unexpected and dangerous!
    pt::Wanderer wanderer(authority, std::move(partition), tree);
  }

  // FIXME: because the first tree put onto the wanderer.trees_ stack
  //        is tree, wanderer will try to destroy it in its destructor
  //        and we're going to get a double free here. see comment in
  //        Wanderer::~Wanderer().

  pll_utree_destroy(tree);
}
