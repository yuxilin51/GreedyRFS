# Author: Xilin Yu @ Apr 17 2019
# Given two input trees T1 and T2 with overlapping leaf sets (in newick form), finds the supertree T* on all leaves that maximizes the bipartitions shared by T* and T1,T2.

import networkx as nx
import networkx.algorithms as alg
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

############### Helper functions ###############


"""
Computes the set of all leaves of T
"""
def leafSet(T):
	return set([v for v,d in T.degree() if d == 1])

"""
Computes the bipartition induced by the edge e in the tree T
"""
def bipartition(T, e):
	# make a copy of T and remove e, dfs on both end-vertex of e to find all leaves in each component  
	T = T.copy()
	T.remove_edge(*e)
	A = set(dfs.dfs_preorder_nodes(T,e[0]))
	B = set(dfs.dfs_preorder_nodes(T,e[1]))
	# returns the bipartition if both sides are non-empty or None otherwise
	return (A,B) if bipartition_is_valid((A,B)) else None

"""
Computes the bipartition restricted to the given subset of leaves X
"""
def restrict_bipartition(pi, X):
	if pi is None:
		return None
	A = pi[0].intersection(X)
	B = pi[1].intersection(X)
	# returns the bipartition if both sides are non-empty or None otherwise
	return (A,B) if bipartition_is_valid((A,B)) else None

def bipartition_is_valid(pi):
	return len(pi[0]) != 0 and len(pi[1]) != 0 

"""
Decides whether the bipartition pi is trivial
"""
def bipartition_is_trivial(pi):
	return len(pi[0]) == 1 or len(pi[1]) == 1


""" 
Returns the set of bipartitions of the tree T restricted to leaves in X
"""
def bipartitions(T, X = None):
	# if X is not given, set it to the whole leaf set 
	if X is None:
		X = leafSet(T)
	# iterate through all edges and add the bipartition induced by that edge if the bipartition is not None
	bipartitions = []
	for e in T.edges():
		pi = restrict_bipartition(bipartition(T,e), X)
		if pi is not None:
			bipartitions.append(pi)
	return bipartitions



""" 
Returns the set of non-trivial bipartitions of the tree T restricted to leaves in X
"""
def non_trivial_bipartitions(T, X = None):
	return [pi for pi in bipartitions(T,X) if not bipartition_is_trivial(pi)]



"""
Returns the (component) subtree in T-e containing the vertex v
"""
def subtree_off_edge(T, e, v):
	T = T.copy()
	T.remove_edge(*e)
	C = alg.components.connected_components(T)
	for c in C:
		if v in c:
			return T.subgraph(c)



"""
Returns the extra subtrees attached to the edge that induces pi in T|_X
"""
def extra_subtrees(T, pi, X = None):
	
	# if X is not given, X is calculated as union of two sides of pi
	if X is None:
		X = pi[0].union(pi[1])

	all_bipartitions = bipartitions(T)
	extra_subtrees = []


def same_bipartition(pi1,pi2):
	if pi1 is None and pi2 is None:
		return True
	elif (pi1 is None) != (pi2 is None):
		return False 
	else:
		return (pi1[0]==pi2[0] and pi1[1] == pi2[1]) or (pi1[0] == pi2[1] and pi1[1] == pi2[0])

def edges_of_bipartition(T, pi):
	for e in T.edges():
		new_pi = restrict_bipartition(bipartition(T,e),X)
		 
	return 

"""

"""
def refine():
	pass



def main():
	T = nx.Graph()
	T.add_nodes_from(['a','b','c','d','e','f','g','ab','abc','dg','dgef','ef'])
	T.add_edges_from([('a','ab'),('b','ab'),('c','abc'),('ab','abc'),('e','ef'),('f','ef'),('ef','dgef'),('d','dg'),('g','dg'),('dg','dgef'),('abc','dgef')])
	# extra_trees = extra_subtrees(T1, ({'a','c'},{'e','f'}))
	b2 = bipartitions(T)
	b1 = non_trivial_bipartitions(T,{'a','d','e','f'})
	b3 = restrict_bipartition(bipartition(T,('ab','abc')),{'a','b','d','e'})
	b4 = restrict_bipartition(bipartition(T,('abc','dgef')),{'a','b','d','e'})
	print(same_bipartition(b3,b4))
	# print(b1)
	# print(b2)
	nx.draw(T, with_labels = True)
	plt.show()

if __name__ == '__main__':
	main()