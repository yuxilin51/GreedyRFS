#!/Users/lethien_96/anaconda/bin/python
"""
Modified by Thien Le in 2019 May for HMMDecomposition
"""
"""
Build subsets using centroid decomposition from PASTA

Written by EKM (molloy.erin.k@gmail.com) in Spring 2018.
"""
import dendropy
import argparse
import copy
from pasta.pastaalignerjob import bisect_tree
from pasta.tree import PhylogeneticTree
import sys

"""
Copied from PASTA
"""

def main(args):
    # Step 1: Decompose tree
    tree = dendropy.Tree.get(path=args.input_tree_file,
                             schema="newick")
    tree.resolve_polytomies(limit=2,
                            update_bipartitions=True)
    
    phy_tree = PhylogeneticTree(tree)

    X = []

    centroid_edge = phy_tree.get_centroid_edge()

    dd_tree = copy.deepcopy(tree)
        

    max_e_w = 0.0
    for e in dd_tree.postorder_edge_iter():
        if e.length == None:
            e.length = 0
        else:
            max_e_w = max(max_e_w, e.length)
        
    max_e_w = 200 * max_e_w

    centroid_edge.length = max_e_w
    a = centroid_edge.head_node
    b = centroid_edge.tail_node
    a_arr = []
    b_arr = [] 
    dd_tree.reroot_at_node(a)
    for leaf in dd_tree.leaf_nodes():
        a_arr.append((leaf.taxon.label, leaf.distance_from_root()))
    a_arr = sorted(a_arr, key = lambda x: x[1])

    for i in range(25):
        X.append(a_arr[i][0])

    dd_tree.reroot_at_node(b)
    for leaf in dd_tree.leaf_nodes():
        b_arr.append((leaf.taxon.label, leaf.distance_from_root()))
    b_arr = sorted(b_arr, key = lambda x: x[1])

    for i in range(25):
        X.append(b_arr[i][0])

    with open(args.output + "/X.lab", "w") as f:
        f.write("\n".join(X))
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="")

    parser.add_argument("-t", "--input_tree_file", type=str,
                        required=True,
                        help="Input tree file")

    parser.add_argument("-o", "--output", type=str,
                        required=True,
                        help="Output prefix")

    main(parser.parse_args())
