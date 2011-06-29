import os,sys
import dendropy
import csv

"""
Quick script for renaming tree labels based on an annotation file
which must be tab-delimited and have the fields
 #SampleID   SampleName

#SampleID will be replaced by SampleName in the tree.
"""

def annotate(tree_filename, annt_filename, output_filename):
	t = dendropy.Tree.get_from_path(tree_filename, 'newick')

	for d in csv.DictReader(open(annt_filename), delimiter='\t'):
		n = t.find_node_with_taxon_label(d["#SampleID"])
		n.taxon.label = d["SampleName"]
		print >> sys.stderr, "renaming {0} to {1}....".format(d["#SampleID"], d["SampleName"])

	with open(output_filename, 'w') as f:
		f.write(t.as_newick_string() + '\n')

if __name__ == "__main__":
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option("-t", "--tree", dest="tree_filename", help="Input tree filename")
	parser.add_option("-a", "--annotation", dest="annt_filename", help="Annotation filename")
	parser.add_option("-o", "--output", dest="output_filename", help="Output tree filename")
	options, args = parser.parse_args()

	annotate(options.tree_filename, options.annt_filename, options.output_filename)
