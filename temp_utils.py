import os
import sys
import numpy as np

def write_DoC(df_list, f, mask=None, delimiter=','):
	"""
	Given <df_list>, a list of DF objects,
	write out a DoC(Depth-of-Coverage) output using get_compressed_vec
	"""
	for df in df_list:
		f.write(df.name + delimiter)
		if mask is None:
			df.get_compressed_vec().tofile(f, sep=delimiter)
		else:
			df.get_compressed_vec()[mask].tofile(f, sep=delimiter)
		f.write('\n')

def combine_df(df_list, new_name=''):
	"""
	Combine all DFs in <df_list> 
	"""
	result = DF(name=new_name, len=df_list[0].len)
	for df in df_list:
		result = result + df
	return result

def create_threshold_mask_for_df_list(df_list, threshold):
	tmp = df_list[0].get_compressed_vec()
	nzs = filter(lambda i: tmp[i] >= threshold, tmp.nonzero()[0])
	for df in df_list[1:]:
		tmp = df.get_compressed_vec()
		nzs = np.lib.arraysetops.intersect1d(nzs, filter(lambda i: tmp[i] >= threshold, tmp.nonzero()[0]))
	return nzs

def annotate_tree_file(tree_filename, annotation_filename, output_filename):
	"""
	Given a tree file and an annotation file (tab-delimited, must have "#SampleID" and "Name" field)
	create a new tree file <output_filename> subbing the sample names with the "Name" field
	"""
	from newick import parse_tree
	tree = parse_tree(open(tree_filename).read())

	import csv
	os.system("cp {0} {1}".format(tree_filename, output_filename))
	sub_dict = dict((r["#SampleID"],r["Name"]) for r in csv.DictReader(open(annotation_filename), delimiter='\t'))

	for le in tree.leaves:
		if le.identifier in sub_dict:
			print >> sys.stderr, "subbing {0} with {1}".format(le.identifier, sub_dict[le.identifier])
			le.identifier = sub_dict[le.identifier]
		else:
			print >> sys.stderr, "ignore {0} not in dict!".format(le.identifier)
		
	with open(output_filename, 'w') as f:
		f.write(str(tree) + '\n')

if __name__ == "__main__":
	annotate_tree_file(*sys.argv[1:])
	sys.exit(-1)
	from DF import DFReader
	import SILVA
	df_filename = sys.argv[1]
	output_filename = sys.argv[2]
	df_list = [df for df in DFReader(open(df_filename))]
	with open(output_filename, 'w') as f:
		write_DoC(df_list, f, mask=SILVA.Ecoli1542_SILVA100)
