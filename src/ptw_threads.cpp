#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include <libpll/pll.h>
#include <model.hpp>
#include <pll_partition.hpp>
#include <pll_util.hpp>

#include "compressed_tree.hpp"
#include "guru.hpp"
#include "options.hpp"
#include "wanderer.hpp"

void WriteTreeTable(pt::TreeTable& trees, const std::string& path)
{
  std::ofstream ofs(path);

  auto table = trees.lock_table();
  for (auto& item : table) {
    ofs << item.second << "\t" << item.first << "\n";
  }
}

int main(int argc, const char* argv[])
{
  pt::Options options = pt::ParseArguments(argc, argv);

  //
  // load tree, sequences, and model info
  //

  std::vector<pll_utree_t*> trees =
      pt::pll::ParseMultiNewick(options.newick_path.c_str());

  if (trees.empty()) {
    throw std::invalid_argument("No valid trees found in " +
                                options.newick_path);
  }

  std::vector<std::string> labels;
  std::vector<std::string> sequences;
  pt::pll::ParseFasta(options.fasta_path, trees[0]->tip_count, labels, sequences);

  pt::CompressedTree::BuildDictionary(labels);

  pt::pll::Model model = pt::pll::ParseRaxmlInfo(options.raxml_path,
                                                 options.rate_categories);

  //
  // initialize the guru
  //

  pt::Guru guru(options, trees, model, labels, sequences);

  //
  // go!
  //

  guru.Start();
  guru.Wait();

  //
  // write output
  //

  if (!options.skip_filtering) {
    // ensure that the good tree table is filtered based on the final threshold
    guru.FilterGoodTreeTable();
  }

  WriteTreeTable(guru.GetGoodTreeTable(), options.good_trees_path);

  if (!options.visited_trees_path.empty()) {
    WriteTreeTable(guru.GetVisitedTreeTable(), options.visited_trees_path);
  }

  if (!options.tested_trees_path.empty()) {
    WriteTreeTable(guru.GetTestedTreeTable(), options.tested_trees_path);
  }

  //
  // clean up and return
  //

  for (auto tree : trees) {
    pll_utree_destroy(tree, pt::pll::cb_erase_data);
  }

  return 0;
}
