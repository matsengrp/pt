[![wercker status](https://app.wercker.com/status/d4b02863ff26dc75109c78a85225cd61/s/master "wercker status")](https://app.wercker.com/project/byKey/d4b02863ff26dc75109c78a85225cd61)

# phylogenetic topographer


## Compiling

Most of `pt`'s dependencies are included as git submodules, so when cloning the repository is it necessary to use the `--recursive` flag, like so:

``` shell
git clone --recursive https://github.com/matsengrp/pt.git
```

A C++11 compiler, CMake 3.0+, and GSL (version 1 or 2) are required for compilation.
To compile, simply run `make` in the project directory.
To run `pt`'s unit tests, run `make test`.
The unit tests usually take between 5-10 minutes to run.

The compiled binary can be found as `_build/src/ptw_threads`.
The binary is statically linked, so no library installation is required.
Simply copy the binary to somewhere in your `PATH` and you're ready to go!
Help for command-line invocation and the available options is available through `ptw_threads --help`.


## Example

This example requires that you have RAxML installed.
We used version 8.2.9.
`cd` into the `test/test-data/hohna_datasets_fasta` directory.

The first step is to run RAxML on the sequence alignment to find a starting point for the `pt` run.
In this case we'll just assume a GTR model with four discrete Gamma rate categories, and leave the rest of the parameters at their defaults.

``` shell
raxmlHPC-AVX -m GTRGAMMA -s DS1.fasta -n DS1
```

(Substitute `raxmlHPC-AVX` with whichever RAxML binary you have available.)
This will generate two files we'll need to run `pt`: `RAxML_info.DS1`, which contains the ML model information, and `RAxML_bestTree.DS1`, which contains the ML tree.

The following command will launch a `pt` run starting at the RAxML ML tree and model, searching out to an offset of -2 units from the log-likelihood of the ML tree.
We set the optimization radius to zero, indicating that on each NNI move, only the branch across which the move is made should be optimized in order to determine if a step should be taken to the resulting tree.

``` shell
ptw_threads --lnl-offset -2 --optimization-radius 0 --raxml-info RAxML_info.DS1 RAxML_bestTree.DS1 DS1.fasta good_trees.tsv
```

The output of the `pt` run, `good_trees.tsv`, is a tab-separated text file where the first column is the log-likelihood of the fully-optimized tree given in the second column.
