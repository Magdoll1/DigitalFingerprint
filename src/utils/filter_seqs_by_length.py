import os
import sys
from Bio import SeqIO

def filter_seqs_by_length(input_filename, output_filename, length_cutoff):
	"""
	Exclude all seqs in the <output> from <input> that have seq length < cutoff
	"""
	excluded = 0
	with open(output_filename, 'w') as f:
		for r in SeqIO.parse(open(input_filename), 'fasta'):
			if len(r.seq) >= length_cutoff:
				f.write(">{id}\n{seq}\n".format(id=r.id, seq=r.seq))
			else:
				excluded += 1
	print("{0} seqs excluded due to short length".format(excluded))

if __name__ == "__main__":
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option("-i", "--input", metavar="FILE")
	parser.add_option("-o", "--output", metavar="FILE")
	parser.add_option("-l", "--length")

	options, args = parser.parse_args()

	filter_seqs_by_length(options.input, options.output, int(options.length))


