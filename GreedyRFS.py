import dendropy as dd
import networkx as nx
import matplotlib.pyplot as plt
import max_bisup_supertree as supertree
import argparse
import os
import sys
import time

"""
Turn the given list of dendropy trees into networkx trees and run max_bisup_tree_many on the list
Turn the resulting networkx supertree back to dendropy tree and write it to file in outdir
"""
def run_max_bisup_suptertree(dd_trees, outfile):
    trees = []
    i = 1
    node_label = dict()
    # node_counter = 0
    for dd_t in dd_trees:
        t = nx.Graph()
        # add nodes and give non-leaves names with number
        for node in dd_t.nodes():
            if node.taxon is not None:
                t.add_node(node.taxon.label)
            else:
                t.add_node('x'+str(i))
                node_label[node] = 'x'+str(i)
                i += 1
        # add edges
        for node in dd_t.nodes():
            tail = node.taxon.label if node.taxon is not None else node_label[node]
            for e in node.child_edge_iter():
                head_taxon = e.head_node.taxon
                head = head_taxon.label if head_taxon is not None else node_label[e.head_node]
                t.add_edge(tail,head)
        # suppress any degree two node (exists at least one b/c of root of dendropy tree)
        suppress_degree_two(t)
        # nx.draw(t, with_labels = True)
        # plt.show()
        # print(i)
        i = supertree.arbitrary_refine(t, i)
        trees.append(t)
    # print("Checking all input" + str(len(dd_trees)))
    # time.sleep(5)
    output_dd_tree = nx_tree_to_dd_tree(supertree.max_bisup_tree_many(trees, i + 1))
    output_dd_tree.write(path = outfile, schema = "newick", suppress_leaf_taxon_labels = True, suppress_leaf_node_labels = False, suppress_internal_node_labels = True)
    # nx_tree_to_dd_tree(trees[0]).write(path = "one", schema = "newick", suppress_leaf_taxon_labels = True, suppress_leaf_node_labels = False, suppress_internal_node_labels = True)
    # nx_tree_to_dd_tree(trees[1]).write(path = "two", schema = "newick", suppress_leaf_taxon_labels = True, suppress_leaf_node_labels = False, suppress_internal_node_labels = True)

"""
Turns a networkx tree into a dendropy tree
"""
def nx_tree_to_dd_tree(T):
    # build dict from string label to dendropy node

    T.add_node('my_root')

    random_node = next(iter(T.nodes()))
    near = next(iter(T.neighbors(random_node)))

    T.remove_edge(random_node, near)
    T.add_edge(random_node, 'my_root')
    T.add_edge(near, 'my_root')

    label_node = dict()
    for v in T.nodes():
        dd_node = dd.Node(label = v)
        label_node[v] = dd_node
    
    # root tree at random node
    # root = next(iter(T.nodes()))
    dd_tree = dd.Tree()
    dd_tree.seed_node = label_node['my_root']
    # print(root)

    # add the edges in the tree
    for v,successors in nx.bfs_successors(T, 'my_root'):
        dd_node = label_node[v]
        for s in successors:
            dd_child = label_node[s]
            dd_node.add_child(dd_child)

    # nx.draw(T, with_labels = True)
    # plt.show()
    return dd_tree

        

"""
Delete degree two nodes and connect its neighbors
"""
def suppress_degree_two(T):
    degree_two_nodes = [v for v,d in T.degree() if d == 2]
    for v in degree_two_nodes:
        nbs = list(T.neighbors(v))
        T.add_edge(nbs[0],nbs[1])
        T.remove_node(v)



"""
Get list newick strings from files and turn into dendropy trees
"""
def main(args):
    with open(args.treesfile,"r") as file:
        dd_trees = [dd.Tree.get(data = l, schema = "newick") for l in file]
    run_max_bisup_suptertree(dd_trees, args.outfile)




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--treesfile", type = str, help = "Path to the file containing all input trees", required = True)
    parser.add_argument("-o", "--outfile", type = str, help = "Path to output file", required = True)
    main(parser.parse_args())
