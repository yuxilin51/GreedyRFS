import dendropy as dd
import networkx as nx
import matplotlib.pyplot as plt
import max_bisup_supertree as supertree
import argparse
import os
import sys

"""
Turn the given list of dendropy trees into networkx trees and run max_bisup_tree_many on the list
Turn the resulting networkx supertree back to dendropy tree and write it to file in outdir
"""
def run_max_bisup_suptertree(dd_trees, outdir):
    trees = []
    i = 1
    node_label = dict()
    for dd_t in dd_trees:
        t = nx.Graph()
        # add nodes and give non-leaves names with number
        for node in dd_t.nodes():
            if node.taxon is not None:
                t.add_node(node.taxon)
            else:
                t.add_node(str(i))
                node_label[node] = str(i)
                i += 1
        # add edges
        for node in dd_t.nodes():
            tail = node.taxon if node.taxon is not None else node_label[node]
            for e in node.child_edge_iter():
                head = e.head_node.taxon
                if head is None:
                    head = node_label[e.head_node]
                t.add_edge(tail,head)
        # suppress any degree two node (exists at least one b/c of root of dendropy tree)
        suppress_degree_two(t)
        # nx.draw(t, with_labels = True)
        # plt.show()
        trees.append(t)
    output_dd_tree = nx_tree_to_dd_tree(supertree.max_bisup_tree_many(trees))
    output_dd_tree.write(path = outdir, schema = "newick")

"""
Turns a networkx tree into a dendropy tree
"""
def nx_tree_to_dd_tree(T):
    return dd.Tree()

"""
Delete degree two nodes and connect its neighbors
"""
def suppress_degree_two(T):
    degree_two_nodes = [v for v,d in T.degree() if d == 2]
    for v in degree_two_nodes:
        nbs = list(T.neighbors(v))
        T.add_edge(nbs[0],nbs[1])
        T.remove_node(v)


def main(args):
    s1 = "(((A,B),C),(D,E));"
    t1 = dd.Tree.get(data = s1, schema = "newick")
    s2 = "((((A,B),C),((D,E),F)),G);"
    t2 = dd.Tree.get(data = s2, schema = "newick")
    run_max_bisup_suptertree([t1,t2],"/home/yuxilin/max_bisup_supertree/output.nwk")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    main(parser.parse_args())