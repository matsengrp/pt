#include "options.hpp"

#include <cstdlib>
#include <memory>
#include <string>

#include <tclap/CmdLine.h>

#include "move_tester/always.hpp"
#include "move_tester/branch_neighborhood_optimizer.hpp"
#include "move_tester/single_branch_optimizer.hpp"

namespace pt {

void ValidateOptions(const Options& options)
{
  if (options.lnl_offset >= 0.0) {
    throw TCLAP::CmdLineParseException(
        "log-likelihood offset must be less than 0",
        "lnl-offset");
  }

  if (options.thread_count < 1) {
    throw TCLAP::CmdLineParseException(
        "number of threads must be at least 1",
        "threads");
  }

  if (options.rate_categories < 1) {
    throw TCLAP::CmdLineParseException(
        "number of rate categories must be at least 1",
        "rate-categories");
  }
}

Options ParseArguments(int argc, const char* argv[])
{
  Options options;

  try {
    TCLAP::CmdLine cmd("pt: the phylogenetic topographer", ' ', "0.1");

    //
    // required value arguments
    //

    TCLAP::ValueArg<double> lnl_offset(
        "",
        "lnl-offset",
        "Negative log-likelihood offset.",
        true,
        0.0,
        "value");
    cmd.add(lnl_offset);

    // TODO: eventually we might want different ways of specifying
    // model parameters, so we could potentially use cmd.xorAdd() with
    // those additional arguments for specifying that they're
    // exclusive
    TCLAP::ValueArg<std::string> raxml_path(
        "",
        "raxml-info",
        "RAxML info file for model parameters.",
        true,
        "",
        "file");
    cmd.add(raxml_path);

    //
    // optional value arguments
    //

    TCLAP::ValueArg<size_t> thread_count(
        "",
        "threads",
        "Number of search threads.",
        false,
        1,
        "value");
    cmd.add(thread_count);

    TCLAP::ValueArg<unsigned int> rate_categories(
        "",
        "rate-categories",
        "Number of discrete Gamma rate categories.",
        false,
        4,
        "value");
    cmd.add(rate_categories);

    TCLAP::ValueArg<unsigned int> poll_ms(
        "",
        "polling-time",
        "Thread polling time, in milliseconds.",
        false,
        100,
        "value");
    cmd.add(poll_ms);

    TCLAP::ValueArg<std::string> visited_trees_path(
        "",
        "visited-trees",
        "Output file for all visited trees.",
        false,
        "",
        "file");
    cmd.add(visited_trees_path);

    TCLAP::ValueArg<std::string> tested_trees_path(
        "",
        "tested-trees",
        "Output file for all tested trees. Note that starting trees "
          "are always visited but never tested (since there is no move "
          "to a starting tree to test), so starting trees will not "
          "appear in this output. Enabling this option significantly "
          "increases runtime memory usage.",
        false,
        "",
        "file");
    cmd.add(tested_trees_path);

    TCLAP::ValueArg<int> optimization_radius(
        "",
        "optimization-radius",
        "Branch length optimization radius. If this argument is not "
          "set, all branches are optimized to determine if a move "
          "should be accepted. A radius of 0 means that only the "
          "branch across which the move is made is optimized. For "
          "values greater than 0, note that the radius is around one of "
          "the nodes at either end of that branch, so a radius of 2 or "
          "greater is required to ensure the other four branches "
          "connected to that branch are optimized.",
        false,
        -1,
        "value");
    cmd.add(optimization_radius);

    //
    // optional switch arguments
    //

    TCLAP::SwitchArg skip_filtering(
        "",
        "skip-filtering",
        "Skip filtering of good trees by the final log-likelihood "
          "threshold once the search is complete. During the search, "
          "if a tree with a higher log-likelihood than the previous "
          "best tree is found, the threshold for good trees is raised. "
          "Normally, the final set of good trees is filtered to "
          "include only trees with log-likelihoods passing the final "
          "threshold. Setting this flag skips that filtering, so that "
          "any tree that passed the threshold at the time it was "
          "evaluated will be included in the results.");
    cmd.add(skip_filtering);

    TCLAP::SwitchArg optimize_models(
        "",
        "optimize-models",
        "After each move, optimize substitution model parameters "
          "in addition to the tree's branch lengths.");
    cmd.add(optimize_models);

    TCLAP::SwitchArg map_mode(
        "",
        "map",
        "Use MAP instead of ML for optimization, assuming an Exp(10) prior.");
    cmd.add(map_mode);

    TCLAP::SwitchArg marginal_mode(
        "",
        "marginal",
        "Use estimated log-marginal likelihood instead of log-likelihood. "
          "Implies --map.");
    cmd.add(marginal_mode);

    //
    // positional arguments
    //

    TCLAP::UnlabeledValueArg<std::string> newick_path(
        "trees-file",
        "Input file containing starting trees (in Newick format, one tree "
          "per line). If multiple trees are specified, searches are "
          "started from the trees in order of descending log-likelihood.",
        true,
        "",
        "trees_file");
    cmd.add(newick_path);

    TCLAP::UnlabeledValueArg<std::string> fasta_path(
        "alignment-file",
        "Input file containing sequence alignment (in FASTA format).",
        true,
        "",
        "alignment_file");
    cmd.add(fasta_path);

    TCLAP::UnlabeledValueArg<std::string> good_trees_path(
        "output-file",
        "Output file for good trees.",
        true,
        "",
        "output_file");
    cmd.add(good_trees_path);

    //
    // parse the arguments
    //

    cmd.parse(argc, argv);

    //
    // load up the options struct
    //

    options.newick_path = newick_path.getValue();
    options.fasta_path = fasta_path.getValue();
    options.raxml_path = raxml_path.getValue();

    options.lnl_offset = lnl_offset.getValue();
    options.thread_count = thread_count.getValue();
    options.rate_categories = rate_categories.getValue();
    options.poll_ms = poll_ms.getValue();
    options.optimization_radius = optimization_radius.getValue();

    options.skip_filtering = skip_filtering.getValue();
    options.optimize_models = optimize_models.getValue();
    options.marginal_mode = marginal_mode.getValue();

    if (options.marginal_mode) {
      options.map_mode = true;
    } else {
      options.map_mode = map_mode.getValue();
    }

    if (options.optimization_radius < 0) {
      options.move_tester = std::make_shared<move_tester::Always>(options);
    } else if (options.optimization_radius == 0) {
      // As currently written, SingleBranchOptimizer is more efficient
      // than BranchNeighborhoodOptimizer with a radius of 0 would be,
      // as it can test the move in-place rather than requiring that
      // the tree be cloned, and only requires a partial traversal
      // afterward instead of a full traversal.
      options.move_tester = std::make_shared<move_tester::SingleBranchOptimizer>(options);
    } else {
      options.move_tester =
          std::make_shared<move_tester::BranchNeighborhoodOptimizer>(options);
    }

    options.good_trees_path = good_trees_path.getValue();
    options.visited_trees_path = visited_trees_path.getValue();
    options.tested_trees_path = tested_trees_path.getValue();

    if (options.tested_trees_path.empty()) {
      options.track_tested_trees = false;
    }

    //
    // validate options
    //

    ValidateOptions(options);
  } catch (TCLAP::ArgException& e) {
    std::cerr << "error: " << e.error() << " for " << e.argId() << "\n";
    std::exit(1);
  }

  return options;
}

} // namespace pt
