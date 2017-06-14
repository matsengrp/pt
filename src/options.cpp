#include "options.hpp"

#include <cstdlib>
#include <string>

#include <tclap/CmdLine.h>

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

  if (options.optimization_radius > 0) {
    throw TCLAP::CmdLineParseException(
        "branch length optimization radius greater than 0 not implemented",
        "optimization-radius");
  }

  // TODO: add checks for input file existence etc?
}

Options ParseArguments(int argc, const char* argv[])
{
  Options options;

  // usage: ptw --lnl-offset <N>
  //            --raxml-info <FILE>
  //            [--threads <N>]
  //            [--visited-trees <FILE>]
  //            <TREE_FILE>
  //            <ALIGNMENT_FILE>
  //            <GOOD_TREES_FILE>

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

    TCLAP::ValueArg<std::string> visited_trees_path(
        "",
        "visited-trees",
        "Output file for all visited trees.",
        false,
        "",
        "file");
    cmd.add(visited_trees_path);

    TCLAP::ValueArg<int> optimization_radius(
        "",
        "optimization-radius",
        "Branch length optimization radius. If this argument is not "
          "set, all branches are optimized to determine if a move "
          "should be accepted. At present the only radius implemented "
          "is 0, meaning only the branch across which the move is made "
          "is optimized.",
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

    //
    // positional arguments
    //

    TCLAP::UnlabeledValueArg<std::string> newick_path(
        "tree-file",
        "Input file containing starting tree (in Newick format).",
        true,
        "",
        "tree_file");
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

    // TODO: this is inelegant, but eventually the try_all_moves
    //       option should be removed when the move testing strategy
    //       is factored out of the wanderer class and the
    //       optimization radius is fully implemented with pll-modules
    options.optimization_radius = optimization_radius.getValue();
    if (options.optimization_radius < 0) {
      options.try_all_moves = true;
    } else {
      options.try_all_moves = false;
    }

    options.skip_filtering = skip_filtering.getValue();

    options.good_trees_path = good_trees_path.getValue();
    options.visited_trees_path = visited_trees_path.getValue();

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
