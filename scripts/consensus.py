import math
import sys
import operator as op
import itertools as it
from collections import defaultdict
from dendropy import Annotation, Tree

def splitPosteriors(trees):
    splits = defaultdict(lambda: 0.0)
    taxa = None
    Z = 0
    for n, (weight, tree) in enumerate(trees):
        expw = math.exp(weight)
        Z += expw
        for split in tree.encode_bipartitions():
            splits[split] += expw
    splits.default_factory = lambda: -math.inf
    for split in splits:
        splits[split] = math.log(splits[split])
    if Z != 1:
        for split in splits:
            splits[split] -= math.log(Z)
    return splits

def consensus(splits, taxon_namespace, p=0.5):
    sorted_splits = sorted(splits.items(), key=op.itemgetter(1), reverse=True)
    majority_splits = map(op.itemgetter(0), it.takewhile(lambda x: x[1] >= math.log(p), sorted_splits))
    t = Tree.from_bipartition_encoding(majority_splits, taxon_namespace)
    for e in t.edges():
        e.annotations.add(Annotation('post', math.exp(splits[e.bipartition])))
    return t

likes = []
trees = []
ll, t = sys.stdin.readline().strip().split()
likes.append(float(ll))
trees.append(Tree.get(data=t, schema='newick'))
ns = trees[0].taxon_namespace
for l in sys.stdin:
    ll, t = l.strip().split()
    likes.append(float(ll))
    trees.append(Tree.get(data=t, schema='newick', taxon_namespace=ns))

ml = max(likes)
trees = zip(map(lambda x: x - ml, likes), trees)
splits = splitPosteriors(trees)
cons = consensus(splits, ns)
cons.write(file=sys.stdout, schema='nexus')
