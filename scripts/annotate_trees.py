import os,sys
import dendropy
import csv
from DigitalFingerprint import DF
"""
Quick script for renaming tree labels based on an annotation file
which must be tab-delimited and have the fields
 #SampleID   SampleName

#SampleID will be replaced by SampleName in the tree.
"""

def annotate_tree(tree_filename, annt_filename, output_filename):
	t = dendropy.Tree.get_from_path(tree_filename, 'newick')

	for d in csv.DictReader(open(annt_filename), delimiter='\t'):
		n = t.find_node_with_taxon_label(d["#SampleID"])
		n.taxon.label = d["SampleName"]
		print >> sys.stderr, "renaming {0} to {1}....".format(d["#SampleID"], d["SampleName"])

	with open(output_filename, 'w') as f:
		f.write(t.as_newick_string() + '\n')

def annotate_df(df_filename, annt_filename, output_filename):
	mapping = {}
	for d in csv.DictReader(open(annt_filename), delimiter='\t'):
		mapping[d["#SampleID"]] = d["SampleName"]

	f = open(output_filename, 'w')
	df_list = [df for df in DF.DFReader(open(df_filename))]
	for df in df_list:
		if df.name in mapping:
			df.name = mapping[df.name]
	w = DF.DFWriter(f)
	w.writes(df_list)
	f.close()

if __name__ == "__main__":
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option("-t", "--tree", dest="tree_filename", help="Input tree filename (dun use -d)")
	parser.add_option("-d", "--df", dest="df_filename", help="Input DF filename (dun use -t)")
	parser.add_option("-a", "--annotation", dest="annt_filename", help="Annotation filename")
	parser.add_option("-o", "--output", dest="output_filename", help="Output tree/di filename")
	options, args = parser.parse_args()

	if options.tree_filename is not None:
		annotate_tree(options.tree_filename, options.annt_filename, options.output_filename)
	elif options.df_filename is not None:
		annotate_df(options.df_filename, options.annt_filename, options.output_filename)
