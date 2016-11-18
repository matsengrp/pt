# Testing pt log-likelihoods vs. RAxML

The output from `Partition::PrintTables()` looks like this:

```
Log Likelihood for (Ref.A1.AU.03.PS1044_Day0.DQ676872,((Ref.A1.RW.92.92RW008.AB253421,Ref.A1.UG.92.92UG037.AB253429),Ref.A2.CD.97.97CDKTB48.AF286238),Ref.A2.CM.01.01CM_1445MV.GU201516); : -3816.22
...
```

The `[multithreading]` test case uses the data in this directory, `test/test-data/five`.
The good trees from the `PrintTables()` output for this test were extracted via copy and paste into a new file as follows:

``` shell
$ cut -d ' ' -f 4 - > pt_good_trees.five
<paste>
```

pt uses the same model parameters for evaluating each tree, so for the purposes of comparison we'd like to have RAxML do the same.
One of RAxML's modes of operation (`-f n`) allows for the evaluation of multiple trees while only estimating model parameters for the first tree in the file.
We make sure the ML tree is the first tree, in hopes that it will arrive at the same parameter values pt is using:

``` shell
$ cat RAxML_bestTree.five pt_good_trees.five > test_trees.five
```

We then run RAxML in `-f n` mode on the trees:

``` shell
$ raxmlHPC-PTHREADS -f n -z test_trees.five -m GTRGAMMA -s five.fasta -n test -p 1
```

Next we need to associate the RAxML log-likelihoods with their respective trees (separated by tabs), and remove the now-unnecessary first line:

``` shell
$ grep Tree-Length RAxML_info.test | cut -d ' ' -f 4 | paste test_trees.five - | tail -n+2 > good_trees.five.raxml

```
