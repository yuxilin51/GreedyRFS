
"""
Modified by Thien Le in 2019 May for HMMDecomposition
"""
"""
Build subsets using centroid decomposition from PASTA

Written by EKM (molloy.erin.k@gmail.com) in Spring 2018.
"""
import dendropy
import argparse
from pasta.pastaalignerjob import bisect_tree
from pasta.tree import PhylogeneticTree
import sys

"""
Copied from PASTA
"""
def bisect_tree(tree):
    e = tree.get_breaking_edge("centroid")
    snl = tree.n_leaves
    tree1, tree2 = tree.bipartition_by_edge(e)
    assert snl == tree1.n_leaves + tree2.n_leaves
    return tree1, tree2

def main(args):
    # Step 1: Decompose tree
    tree = dendropy.Tree.get(path=args.input_tree_file,
                             schema="newick")
    tree.resolve_polytomies(limit=2,
                            update_bipartitions=True)
    tree = PhylogeneticTree(tree)
    t1, t2 = bisect_tree(tree)
    trees = [t1, t2]

    # Step 2: Write out leaf subsets
    # i = 1
    i = 0

    keep1 = t1.leaf_node_names()
    with open(args.output + "/A.lab", "w") as f:
        f.write("\n".join(keep1))
    
    keep2 = t2.leaf_node_names()
    with open(args.output + "/B.lab", "w") as f:
        f.write("\n".join(keep2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-t", "--input_tree_file", type=str,
                        required=True,
                        help="Input tree file")

    parser.add_argument("-o", "--output", type=str,
                        required=True,
                        help="Output prefix")

    main(parser.parse_args())
