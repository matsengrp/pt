#!/usr/bin/env python3
"""
Usage: consensus.py < weighted.trees > consensus.nex
Reads a stream of white-space delimited log-likelihoods and newick trees from stdin.
Prints the credibility-annotated majority-rule consensus tree in NEXUS format to stdout.
"""

import math
import sys
import operator as op
import itertools as it
from collections import defaultdict
from dendropy import Annotation, Tree

def splitCredibilities(weighted_trees):
    """Takes an iterable of (weight, tree) tuples and returns a dictionary of splits and their credibilities.
    The weights need not be normalized.
    """

    splits = defaultdict(lambda: 0.0)
    Z = 0 # normalizing constant

    for weight, tree in weighted_trees:
        Z += weight
        for split in tree.encode_bipartitions():
            splits[split] += weight

    # Normalize credibilities
    for split in splits:
        splits[split] /= Z

    return splits

def consensus(splits, taxon_namespace, q=0.5):
    """Creates a consensus tree from a dictionary of splits w/ credibilities (default is majority-rule consensus).
    Splits are added to the tree in order of decreasing credibility. A split is added if:
    1. Its credibility is >= q
    2. It is compatible with all preceding splits
    The edges of the returned tree are annotated with their split's credibilities.
    """

    # Order and select splits
    sorted_splits = sorted(splits.items(), key=op.itemgetter(1), reverse=True)
    selected_splits = map(op.itemgetter(0), it.takewhile(lambda x: x[1] >= q, sorted_splits))

    # Dendropy handles dropping of incompatible splits
    tree = Tree.from_bipartition_encoding(selected_splits, taxon_namespace)

    # Annotate edges with credibilities
    for edge in tree.edges():
        edge.annotations.add(Annotation('cred', splits[edge.bipartition]))

    return tree

def parse_line(l, taxon_namespace=None):
    l = l.strip().split()
    # First whitespace delimits weight from tree
    return float(l[0]), Tree.get(data=''.join(l[1:]), schema='newick', taxon_namespace=taxon_namespace)

likes = []
trees = []

# Read the first tree to establish the namespace
ll, t = parse_line(sys.stdin.readline())
likes.append(ll)
trees.append(t)
ns = trees[0].taxon_namespace

# Read the rest of the trees
for l in sys.stdin:
    ll, t = parse_line(l, ns)
    likes.append(ll)
    trees.append(t)

# Normalize by largest likelihood to avoid numerical over/under-flow
ml = max(likes)
trees = zip(map(lambda x: math.exp(x - ml), likes), trees)

splits = splitCredibilities(trees)
cons = consensus(splits, ns)
cons.write(file=sys.stdout, schema='nexus')
