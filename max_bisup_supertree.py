# Author: Xilin Yu @ Apr 17 2019
# Given two input trees T1 and T2 with overlapping leaf sets (in newick form), finds the supertree T* on all leaves that maximizes the bipartitions shared by T* and T1,T2.
import dendropy as dd
from random import random
import math
import networkx as nx
import networkx.algorithms as alg
import networkx.algorithms.traversal.depth_first_search as dfs 
import matplotlib.pyplot as plt


"""
Heuristic for computing split fit supertree (and thus RF supertree) of the list of given trees under the given greedy ordering 
"""
def max_bisup_tree_many(trees, node_counter, ordering = "max_shared_leaves"):
	n = len(trees)
	used_tree_idxs = set()
	if ordering is "max_shared_leaves":
		leaf_sets = [leaf_set(t) for t in trees]
		shared_leaves = dict()
		for i in range(n):
			for j in range(i+1, n):
				shared_leaves[(i,j)] = len(leaf_sets[i] & leaf_sets[j])
		for k in range(n-1):
			idx1,idx2 = max(shared_leaves, key = shared_leaves.get)
			new_t, node_counter = max_bisup_tree_two(trees[idx1],trees[idx2], node_counter)
			new_t_dd = nx_tree_to_dd_tree(new_t)
			used_tree_idxs.update({idx1,idx2})
			# compute num shared leaves between newtree and all existing but not used trees and add to dictionary
			trees.append(new_t)
			leaf_sets.append(leaf_set(new_t))
			for l in range(n+k):
				if l in used_tree_idxs:
					continue
				shared_leaves[(l,n+k)] = len(leaf_sets[l] & leaf_sets[n+k])
			# delte all entries related to trees at idx1 and idx2 in dictionary
			for l in range(n+k+1):
				if (idx1,l) in shared_leaves:
					del shared_leaves[(idx1,l)]
				if (l, idx1) in shared_leaves:
					del shared_leaves[(l,idx1)]
				if (idx2,l) in shared_leaves:
					del shared_leaves[(idx2,l)]
				if (l, idx2) in shared_leaves:
					del shared_leaves[(l,idx2)]

	elif ordering is "max_shared_bipartitions":
		pass
	else:
		pass

	return trees[2*n-2]


"""
Computes and returns the maximum bipartition support supertree of the two given trees
"""
def max_bisup_tree_two(T1, T2, node_counter):

	if not nx.is_tree(T1):
		print("Input T1 for max_bisup_tree method is not a tree.")
	elif not alg.is_tree(T2):
		print("Input T2 for max_bisup_tree method is not a tree.")  
	S1 = leaf_set(T1)
	S2 = leaf_set(T2)
	S = S1 | S2
	X = S1 & S2

	CT1X = bipartitions(T1,X)
	CT2X = bipartitions(T2,X)

	# make changes to bipartitions in CT2X such that if it also shows up in CT1X, the order of the tuple follows the one drom CT1X
	# this makes sure that we can directly do set union and intersection on CT1X and CT2X
	# and that the ordered trees respect the same order
	changes = set()
	for pi in CT2X:
		if pi not in CT1X and equiv_bipar(pi) in CT1X:
			changes.add(pi)
	for pi in changes:
		CT2X.remove(pi)
		CT2X.add(equiv_bipar(pi))

	bipar_subtrees = dict()
	weight = dict()
	all_extra_subtrees = set()
	extra_subtrees_to_pi = dict()
	
	# compute ordered_subtrees associated with each bipartition and the weight of each bipartition
	for pi in CT1X:
		bipar_subtrees[pi] = ordered_subtrees(T1,pi,extra_subtrees_to_pi)
		weight[pi] = len(edges_of_bipartition(T1,pi))
		all_extra_subtrees.update(bipar_subtrees[pi])
	for pi in CT2X:
		if pi not in CT1X:
			bipar_subtrees[pi] = ordered_subtrees(T2,pi,extra_subtrees_to_pi)
			weight[pi] = len(edges_of_bipartition(T2,pi))
		else:
			bipar_subtrees[pi].extend(ordered_subtrees(T2,pi,extra_subtrees_to_pi))
			weight[pi] = weight[pi] + len(edges_of_bipartition(T2,pi))			
		all_extra_subtrees.update(bipar_subtrees[pi])
	
	# compute T = T_init, which is a star on leaf set of X with all extra subtrees attached to the center
	T = nx.Graph()
	T.add_node('x'+str(node_counter))
	for x in X:
		T.add_edge('x'+str(node_counter),x)
	for r,t in all_extra_subtrees:
		T.update(t)
		T.add_edge('x'+str(node_counter),r)
	

	# vtx_bipars is a dictionary between a vertice v and a set of bipartitions whose addition requires refinement at v
	vtx_bipars = dict()
	# bipar_vtx is a dictionary between a biparition pi and the vertex that should be refined for the addition of pi
	bipar_vtx = dict()
	# at beginning, vertex '0' associated with every bipartition and every bipartition points to '0'
	vtx_bipars['x'+str(node_counter)] = CT1X | CT2X
	for pi in CT1X | CT2X:
		bipar_vtx[pi] = 'x'+str(node_counter)
	node_counter += 1

	# compute maximum independent set
	G = incompatibility_graph(CT1X, CT2X)
	I = max_ind_set(G, CT1X, CT2X, weight)

	# refine the tree with bipartitions in I and in intersection of CT1X and CT2X
	i = 0
	for pi in I | (CT1X & CT2X):
		node_counter = refine(T, pi, vtx_bipars, bipar_vtx, bipar_subtrees, node_counter, extra_subtrees_to_pi)

	# arbitrarily refine the tree if there is polytomy
	node_counter = arbitrary_refine(T, node_counter)

	return T, node_counter


"""
Constructs and returns the bipartite incompatibility graph of bipartitions in C1 - C2 and C2 - C1, i.e.,  the symmetric difference of C(T1|_X) and C(T2|_X).  
"""
def incompatibility_graph(C1,C2):
	G = nx.DiGraph()
	G.add_nodes_from(C1-C2)
	G.add_nodes_from(C2-C1)
	for pi1 in C1-C2:
		for pi2 in C2-C1:
			if not bipartitions_are_compatible(pi1,pi2):
				G.add_edge(pi1, pi2, capacity = math.inf)
	return G



"""
Computes the maximum independent set in the given directed bipartite graph.
First, add a source s and an edge from s to each of C1-C2, then, add a sink t and an edge from each of C2-C1 to t.
Each new edge has capacity
Find min-cut max flow on resulting graph and return 
"""
def max_ind_set(G,C1,C2, weight):
	G.add_nodes_from(['s','t'])
	for pi in C1-C2:
		G.add_edge('s',pi,capacity = weight[pi])
	for pi in C2-C1:
		G.add_edge(pi,'t',capacity = weight[pi])
	cut_value, cut_partition = nx.minimum_cut(G, 's', 't')
	return (cut_partition[0] & (C1-C2)) | (cut_partition[1] & (C2-C1))


"""

"""
def refine(T, pi, vtx_bipars, bipar_vtx, bipar_subtrees, node_counter, extra_subtrees_to_pi):
	A = pi[0]
	B = pi[1]
	
	# Case 1: pi is trivial, do not split vertex, but attach all subtrees in order
	if bipar_is_trivial(pi):
		# find va and vb 
		if len(A) == 1:
			va = next(iter(A))
			vb = bipar_vtx[pi]
		else:
			va = bipar_vtx[pi]
			vb = next(iter(B))

		# remove edge (va, vb)
		T.remove_edge(va,vb)
		# create a new vertex for each subtree connect through va to vb
		prev = va
		attach_point = bipar_vtx[pi]
		for r,t in bipar_subtrees[pi]:
			T.remove_edge(r, attach_point)
			# current = va+r
			current = 'x'+str(node_counter)
			node_counter += 1
			T.add_node(current)
			T.add_edge(prev,current)
			T.add_edge(current,r)
			prev = current
		# connect last prev to vb
		T.add_edge(prev, vb)		
		
		# delete pi from vtx_bipars and bipar_vtx (may not be necessary?)
		vtx_bipars[attach_point].remove(pi)
		del bipar_vtx[pi]

	# Case 2: pi is not trivial
	else:
		# find vtx to split as original 
		original = bipar_vtx[pi]

		# make two copies of original named va and vb
		va = 'x'+str(node_counter)
		node_counter += 1
		vb = 'x'+str(node_counter)
		node_counter += 1

		# create a new vertex for each subtree associated with pi connect through va to vb
		prev = va
		for r,t in bipar_subtrees[pi]:
			T.remove_edge(r, original)
			# current = va+r
			current = 'x'+str(node_counter)
			node_counter += 1
			T.add_node(current)
			T.add_edge(prev,current)
			T.add_edge(current,r)
			prev = current
		# connect last prev to vb
		T.add_edge(prev, vb)

		vtx_bipars[va] = set()
		vtx_bipars[vb] = set()
		incompatible_bipars = set()
		# update the pointers of all bipartitions associated with original which are compatible with pi
		for bipar in vtx_bipars[original]:
			if not same_bipartition(bipar, pi) and bipartitions_are_compatible(bipar, pi):
				if bipar[0] <= A or bipar[1] <= A:
					vtx_bipars[va].add(bipar)
					bipar_vtx[bipar] = va
				elif bipar[0] <= B or bipar[1] <= B:
					vtx_bipars[vb].add(bipar)
					bipar_vtx[bipar] = vb
				else:
					print("Something is wrong because ", bipar, " is compatible with ", pi, " but its sides are not subsets of A or B")
			elif not bipartitions_are_compatible(bipar,pi):
				incompatible_bipars.add(bipar)

		# partition other subtrees attached to the vertex
		for bipar in vtx_bipars[va]:
			for r,t in bipar_subtrees[bipar]:
				T.remove_edge(r, original)
				T.add_edge(r, va)
		for bipar in vtx_bipars[vb]:
			for r,t in bipar_subtrees[bipar]:
				T.remove_edge(r, original)
				T.add_edge(r, vb)
		for bipar in incompatible_bipars:
			for r,t in bipar_subtrees[bipar]:
				T.remove_edge(r, original)
				T.add_edge(r, va)


		# partition neighbors of original into two sets, 
		#connect those in components containing leaves of A to va and connect those in components containing leaves B to vb
		neighbors = set(T.neighbors(original)).copy()
		T.remove_node(original)
		
		for n in neighbors:
			components = nx.connected_components(T)
			for c in components:
				if n in c:
					# if the component also contains a vertex from A, connect it to va
					if len(c & A) != 0:
						T.add_edge(n, va)
					# else if the component also contains a vertex from B, connect it to vb
					elif len(c & B) != 0:
						T.add_edge(n, vb)
					# ow, n is the root of an extra subtree 
					elif extra_subtrees_to_pi[n][0] <= A or extra_subtrees_to_pi[n][1] <= A:						
						T.add_edge(n, va)
					elif extra_subtrees_to_pi[n][0] <= B or extra_subtrees_to_pi[n][1] <= B:						
						T.add_edge(n, vb)
					else:
						T.add_edge(n, vb)
					break

		# delete original as a key in vtx_bipars
		del vtx_bipars[original]

		# delete pi as a key in bipar_vtx (may not be necessary?)
		del bipar_vtx[pi]

	return node_counter

"""
Refine the given tree at any node with degree > 3 arbitrarily 
(except that both copies of the node need to have at least 2 neighbors other than each other)
"""
def arbitrary_refine(T, node_counter):
	polytomies = {v for v,d in T.degree() if d > 3}
	while len(polytomies) != 0:
		v = next(iter(polytomies))
		va = 'x'+str(node_counter)
		node_counter += 1

		vb = 'x'+str(node_counter)
		node_counter += 1

		T.add_nodes_from([va,vb])
		for n in T.neighbors(v):
			if T.degree(va) < 2:
				T.add_edge(n,va)
			elif T.degree(vb) < 2:
				T.add_edge(n,vb)
			elif random() < 0.5:
				T.add_edge(n,va)
			else:
				T.add_edge(n,vb)
		T.add_edge(va,vb)
		T.remove_node(v)
		polytomies.remove(v)
		if T.degree(va) > 3:
			polytomies.add(va)
		if T.degree(vb) > 3:
			polytomies.add(vb)

	return node_counter


############### Helper functions ###############


"""
Computes the set of all leaves of T.
"""
def leaf_set(T):
	return set([v for v,d in T.degree() if d == 1])




"""
Computes the bipartition induced by the edge e in the tree T.
"""
def bipartition(T, e):
	# make a copy of T and remove e, dfs on both end-vertex of e to find all leaves in each component  
	T = T.copy()
	T.remove_edge(*e)
	A = set(dfs.dfs_preorder_nodes(T,e[0]))
	B = set(dfs.dfs_preorder_nodes(T,e[1]))
	# returns the bipartition if both sides are non-empty or None otherwise
	return (frozenset(A),frozenset(B)) if bipartition_is_valid((A,B)) else None




"""
Computes the bipartition restricted to the given subset of leaves X.
"""
def restrict_bipartition(pi, X):
	if pi is None:
		return None
	A = pi[0].intersection(X)
	B = pi[1].intersection(X)
	# returns the bipartition if both sides are non-empty or None otherwise
	return (frozenset(A),frozenset(B)) if bipartition_is_valid((A,B)) else None



"""
Returns whether the bipartition has both sides non-empty.
"""
def bipartition_is_valid(pi):
	return len(pi[0]) != 0 and len(pi[1]) != 0 




"""
Returns whether the bipartition pi is trivial, i.e., at least one side has only 1 element.
"""
def bipar_is_trivial(pi):
	return len(pi[0]) == 1 or len(pi[1]) == 1



"""
Returns whether two bipartitions of the same leafset are compatible.
Return True if A1 is a subset of A2 of B2 (and consequently B1 is a superset of B2 or A2), 
	or if A2 or B2 is a subset of A1 (and consequently B2 or A2 is a superset of B1).
"""
def bipartitions_are_compatible(pi1, pi2):
	A1 = pi1[0]
	A2 = pi2[0]
	B2 = pi2[1]
	return A1 <= A2 or A2 <= A1 or A1 <= B2 or B2 <= A1



""" 
Returns the set of bipartitions of the tree T restricted to leaves in X.
"""
def bipartitions(T, X = None):
	# if X is not given, set it to the whole leaf set 
	if X is None:
		X = leaf_set(T)
	# iterate through all edges and add the bipartition induced by that edge if the bipartition is not None
	bipartitions = set()
	for e in T.edges():
		pi = restrict_bipartition(bipartition(T,e), X)
		if pi is not None:
			pi_reverse = (pi[1],pi[0])
			if pi not in bipartitions and pi_reverse not in bipartitions:
				bipartitions.add(pi)
	return bipartitions



""" 
Returns the set of non-trivial bipartitions of the tree T restricted to leaves in X.
"""
def non_trivial_bipartitions(T, X = None):
	return set([pi for pi in bipartitions(T,X) if not bipar_is_trivial(pi)])



"""
Returns if two bipartitions are the same.
"""
def same_bipartition(pi1,pi2):
	if pi1 is None and pi2 is None:
		return True
	elif (pi1 is None) != (pi2 is None):
		return False 
	else:
		return (pi1[0]==pi2[0] and pi1[1] == pi2[1]) or (pi1[0] == pi2[1] and pi1[1] == pi2[0])

"""
Returns an equivalent bipartition where the order of two sides in the tuple is swapped
"""
def equiv_bipar(pi):
	return (pi[1],pi[0])

"""
Returns the (component) subtree in T-e containing the vertex v.
"""
def subtree_off_edge(T, e, v):
	T = T.copy()
	T.remove_edge(*e)
	C = alg.components.connected_components(T)

	for c in C:
		if v in c:
			return T.subgraph(c)



"""
Returns an ordered list of extra subtrees on the path of edges which induce pi = (A,B) in T|_X where X = leaves of pi
The first extra subtree is closest to A
"""
def ordered_subtrees(T, pi, extra_subtrees_to_pi):
	#count the appearance of nodes in the set of edges and the ones showing up once are the two ends s and t
	nodes_count = dict()
	for e in edges_of_bipartition(T,pi):
		for v in e:
			if v not in nodes_count:
				nodes_count[v] = 1
			else:
				nodes_count[v] = nodes_count[v] + 1
	ends = [v for v,count in nodes_count.items() if count == 1]
	a = next(iter(pi[0]))
	dist0 = nx.dijkstra_path_length(T,a,ends[0])
	dist1 = nx.dijkstra_path_length(T,a,ends[1])
	if dist0 < dist1:
		s = ends[0]
		t = ends[1]
	else:
		s = ends[1]
		t = ends[0]
	# find the path between s and t, excluding s and t
	inner_nodes = nx.shortest_path(T, s, t)[1:-1]
	# compute edge node pairs of (edge between inner node to the root of the extra subtree, the root of the extra subtree)
	# where roots of the extra subtrees are the neighbors of inner node that do not show up in nodes_on_path 
	edge_node_pairs = []
	for node in inner_nodes:
		for other in T.neighbors(node):
			if other not in nodes_count:
				edge_node_pairs.append(((node,other),other))
	trees = [(v, subtree_off_edge(T,e,v)) for e,v in edge_node_pairs]
	for (v,t) in trees:
		extra_subtrees_to_pi[v] = pi
	return trees


"""
Returns a set of edges which induce bipartitions such that when restricted to the leaf set of pi is equivalent to pi.
"""
def edges_of_bipartition(T, pi):
	X = pi[0].union(pi[1])
	edges = []
	for e in T.edges():
		new_pi = restrict_bipartition(bipartition(T,e),X)
		if same_bipartition(pi, new_pi):
			edges.append(e)
	return edges


"""
Turns a networkx tree into a dendropy tree
"""
def nx_tree_to_dd_tree(T):

    
    # build dict from string label to dendropy node   
    label_node = dict()
    for v in T.nodes():
        dd_node = dd.Node(label = v)
        label_node[v] = dd_node
    
    # root tree at random node
    root = next(iter(T.nodes()))
    dd_tree = dd.Tree()
    dd_tree.seed_node = label_node[root]

    # add the edges in the tree
    for v,successors in nx.bfs_successors(T, root):
        dd_node = label_node[v]
        for s in successors:
            dd_child = label_node[s]
            dd_node.add_child(dd_child)
    return dd_tree


