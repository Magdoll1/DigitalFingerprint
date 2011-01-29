import sys
from sets import ImmutableSet
from newick import parse_tree
from collections import defaultdict

def flatten_help(tree, acc):
	_x = list(tree.leaves_identifiers)
	_x.sort()
	if len(_x) == 1: return 
	acc.append(tuple(_x))
	if len(tree.leaves) <= 2: return
	else:
		flatten_help(tree.get_edges()[0][0], acc)
		flatten_help(tree.get_edges()[1][0], acc)

def flatten(tree):
	acc = []
	flatten_help(tree, acc)
	return acc

def evaluate_jackknifing(real_tree_filename, subsampled_tree_filename):
	real_node_hit_count = {}
	fale = defaultdict(lambda: 0)
	with open(real_tree_filename) as f:
		real_tree = parse_tree(f.read())
		for real_node in flatten(real_tree):
			real_node_hit_count[real_node] = 0
	with open(subsampled_tree_filename) as f:
		for line in f:
			tree = parse_tree(line)
			for node in flatten(tree):
				if node in real_node_hit_count:
					real_node_hit_count[node] += 1
				else:
					fale[node] += 1
	
	return real_node_hit_count, fale


h,f = evaluate_jackknifing(*sys.argv[1:])
