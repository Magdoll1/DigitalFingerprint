import os, sys, glob
from DigitalFingerprint import Read
from DigitalFingerprint.DF import DFWriter
from optparse import OptionParser

from DigitalFingerprint.utils import versatile_open
open = versatile_open.versatile_open

parser = OptionParser()
parser.add_option("-p", "--pattern", dest="pattern", help="pattern for list of input files")
parser.add_option("-o", "--output", dest="output", help="output filename (suffix should be .DF)")
parser.add_option("-r", "--ref_gap_map", dest="ref_gap_map", help="Gap map file of reference alignment (default: Silva/SILVA100_justHumanCrap.1crap_masked.V3region.fna.gap_map.bz2")
parser.add_option("-l", "--ref_aln_len", dest="ref_aln_len", type=int, help="Reference alignment length, must be given if not using default gap map. Default 50,000.")

options, args = parser.parse_args()

pattern = options.pattern
output_df_filename = options.output

if pattern is None or output_df_filename is None:
	print >> sys.stderr, "Must provide input file name or pattern AND output filename!"
	sys.exit(-1)

print >> sys.stderr, "Reading RefMap....this may take a while"
if options.ref_gap_map is None:
	refmap = Read.RefMap('../data/SILVA104.fece_augmented.fasta.gap_map.bz2', aln_length=50000)
else:
	refmap = Read.RefMap(options.ref_gap_map, options.ref_aln_len)

f = open(output_df_filename, 'w')
dfwriter = DFWriter(f)

for file in glob.iglob(pattern):
	# file name ex: O2.UC-1_090112.fq.bowtied
	name = os.path.basename(file)
	print >> sys.stderr, "reading {0} for DF writing....".format(file)
	readdict = Read.ReadsDict(refmap)
	readdict.read_bowtie_output(file)
	readdf = Read.ReadDF(name, refmap)
	for read in readdict:
		readdf.add_read_to_vec(read)
	dfwriter.write(readdf)

f.close()
