# Author: Xilin Yu @ Apr 17 2019
# Given two input trees T1 and T2 with overlapping leaf sets (in newick form), finds the supertree T* on all leaves that maximizes the bipartitions shared by T* and T1,T2.

import networkx as nx
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


	

def refine():
	pass

def main():
	T1 = nx.Graph()
	T1.add_nodes_from(['a','b','c','d','e','f','ab','abc','def','ef'])
	T1.add_edges_from([('a','ab'),('b','ab'),('c','abc'),('ab','abc'),('e','ef'),('f','ef'),('ef','def'),('d','def'),('abc','def')])
	# T1R = nx.subgraph(T1,['a','b','e','f'])
	T1.copy()
	T1.remove_nodes_from(['c','d'])
	nx.draw(T1, with_labels = True)
	plt.show()

if __name__ == '__main__':
	main()
