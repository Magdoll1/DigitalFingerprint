import os, sys
import dendropy
"""
Dendropy's current sum_tree.py seems a bit odd....(it flattened one of my branches)
so I'm doing my simple implementation to double-check
"""

def sum_tree(true_tree_filename, tree_list_filename):
	"""
	vanilla version of sumtree.py from Dendropy
	right now it is slow and stupid...just puts the number of times the subtree
	is seen in the tree_list (not percentage! just abs count!) in the edge
	as a string support:edge_length and return the true tree....
	"""
	#result = []
	true_tree = dendropy.Tree.get_from_path(true_tree_filename, 'newick')
	#taxon_set = true_tree.infer_taxa()
	tree_list = dendropy.TreeList().get_from_path(tree_list_filename, 'newick', taxon_set=taxon_set)
#	split_distribution = dendropy.treesplit.SplitDitribution(taxon_set=taxon_set)
	for p in true_tree.preorder_node_iter():
		if p.is_leaf(): continue
		t_p = dendropy.Tree.get_from_string(p.as_newick_string(), 'newick')
		found = 0
		for t in tree_list:
			for q in t.preorder_node_iter():
				if q.is_leaf(): continue
				t_q = dendropy.Tree.get_from_string(q.as_newick_string(), 'newick')
				if t_p.symmetric_difference(t_q) == 0:
					found += 1
					break
		#result.append((p, found))
		p.edge.length = "{support}:{edge}".format(support=found, edge=p.edge.length)
	return true_tree
