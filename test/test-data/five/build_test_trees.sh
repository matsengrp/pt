#!/bin/bash

set -eu

rm -f RAxML_*.test_trees

../../../_build/src/ptw RAxML_bestTree.five five.fasta RAxML_info.five -82.53 ptw_good.tsv ptw_all.tsv
sort -n -r ptw_good.tsv | cut -f 2 > test_trees.five
raxmlHPC-PTHREADS -f n -z test_trees.five -m GTRGAMMA -s five.fasta -n test_trees -p 1
grep "Tree-Length" RAxML_info.test_trees | cut -d ' ' -f 4 | paste test_trees.five - > good_trees.five.raxml
