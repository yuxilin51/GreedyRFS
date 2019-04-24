# Author: Xilin Yu @ Apr 17 2019
# Given two input trees T1 and T2 with overlapping leaf sets (in newick form), finds the supertree T* on all leaves that maximizes the bipartitions shared by T* and T1,T2.

import dendropy
import networkx as nx
import networkx.algorithms.traversal.depth_first_search as dfs 
import matplotlib.pyplot as plt

def max_bisup_tree(T1, T2):
	if not nx.isTree(T1):
		print("Input T1 for max_bisup_tree method is not a tree.")
	elif not nx.isTree(T2):
		print("Input T2 for max_bisup_tree method is not a tree.")  
	S1 = set(T1.nodes())
	S2 = set(T2,nodes())
	S = S1.union(S2)
	X = S1.intersection(S2)



"""
Computes the bipartition induced by the edge e in the tree T, restricted to the subset X of leaves
"""
def bipartition(T, e, X = None):

	# if no leaf subset is given, compute the set of all leaves as X
	if X is None:
		X = set([v for v,d in T.degree() if d == 1])

	# make a copy of T and remove e, dfs on both end-vertex of e to find all leaves in each component  
	T = T.copy()
	T.remove_edge(*e)
	A = set(dfs.dfs_preorder_nodes(T,e[0])).intersection(X)
	B = set(dfs.dfs_preorder_nodes(T,e[1])).intersection(X)

	# returns the bipartition if both sides are non-empty or None otherwise
	return (A,B) if not len(A) == 0 and not len(B) == 0 else None



"""
Decides whether the bipartition pi is trivial
"""
def trivial_bipartition(pi):
	if len(pi[0]) == 1 or len(pi[1]) == 1:
		return True
	else:
		return False



""" 
Returns the set of bipartitions of the tree T restricted to leaves in X
"""
def bipartitions(T,X  = None):

	# if no leaf subset is given, compute the set of all leaves as X
	if X is None:
		X = set([v for v,d in T.degree() if d == 1])

	# iterate through all edges and add the bipartition induced by that edge if the bipartition is not None
	bipartitions = []
	for e in T.edges():
		pi = bipartition(T,e,X)
		if pi is not None:
			bipartitions.append(pi)
	return bipartitions



""" 
Returns the set of non-trivial bipartitions of the tree T restricted to leaves in X
"""
def non_trivial_bipartitions(T, X = None):

	return [pi for pi in bipartitions(T,X) if not trivial_bipartition(pi)]

"""
Returns the extra subtrees attached to the edge that induces pi in T|_X
"""
def extra_subtrees(T, pi, X = None):
	if X is None:
		X = pi[0].union(pi[1])
	



def edges_of_restricted_bipar(T, pi):

	pass	

"""

"""
def refine():
	pass



def main():
	T1 = nx.Graph()
	T1.add_nodes_from(['a','b','c','d','e','f','ab','abc','def','ef'])
	T1.add_edges_from([('a','ab'),('b','ab'),('c','abc'),('ab','abc'),('e','ef'),('f','ef'),('ef','def'),('d','def'),('abc','def')])
	extra_trees = extra_subtrees(T1, ({'a','c'},{'e','f'}))
	# b1 = non_trivial_bipartitions(T1,{'a','d','e','f'})
	# print(b1)
	# nx.draw(T1, with_labels = True)
	# plt.show()

if __name__ == '__main__':
	main()
