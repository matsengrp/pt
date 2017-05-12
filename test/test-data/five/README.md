# Testing pt log-likelihoods vs. RAxML

This process is automated by the `build_test_trees.sh` script.

`ptw` outputs tab-separated files that look like this:

```
-3816.22	(Ref.A1.AU.03.PS1044_Day0.DQ676872,((Ref.A1.RW.92.92RW008.AB253421,Ref.A1.UG.92.92UG037.AB253429),Ref.A2.CD.97.97CDKTB48.AF286238),Ref.A2.CM.01.01CM_1445MV.GU201516);
...
```

The `[wanderer]` test case uses the data in this directory, `test/test-data/five`.
Trees and log-likelihoods were obtained from `ptw` using the following command in this directory:

``` shell
$ ../../../_build/src/ptw RAxML_bestTree.five five.fasta RAxML_info.five -82.53 ptw_good.tsv ptw_all.tsv
```

pt uses the same model parameters for evaluating each tree, so for the purposes of comparison we'd like to have RAxML do the same.
One of RAxML's modes of operation (`-f n`) allows for the evaluation of multiple trees while only estimating model parameters for the first tree in the file.
The idea is that RAxML will arrive at the same parameter values for this first tree that pt reads from the `RAxML_info.five` file.
We make sure the ML tree is the first tree in the file using a descending numerical sort (since the log-likelihoods are the first field in the `ptw` output files) and then cut out the Newick strings:

``` shell
$ sort -n -r ptw_good.tsv | cut -f 2 > test_trees.five
```

We then run RAxML in `-f n` mode on the trees:

``` shell
$ raxmlHPC-PTHREADS -f n -z test_trees.five -m GTRGAMMA -s five.fasta -n test_trees -p 1
```

Next we need to associate the RAxML log-likelihoods with their respective trees (separated by tabs):

``` shell
$ grep "Tree-Length" RAxML_info.test_trees | cut -d ' ' -f 4 | paste test_trees.five - > good_trees.five.raxml
```
