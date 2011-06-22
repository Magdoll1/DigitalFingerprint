import os, sys, glob
import numpy as np
from DigitalFingerprint.DF import DFReader
from DigitalFingerprint import SILVA
from DigitalFingerprint.clustering import Cluster
from optparse import make_option, OptionParser

def main(args=None):
	ecoli_map = SILVA.Ecoli1542_SILVA100
	
	parser = OptionParser(option_list=[ \
			make_option("-f", "--df-filename", dest="df_filename"), \
			make_option("-r", "--ecoli-range", dest="ecoli_range", help="specify the (1-based) E.coli position range to use ex: -r 200,300"), \
			make_option("-a", "--all-positions", dest="ecoli_only", action="store_false", default=True, help="use all <len> positions, ignore ecoli range"),\
			make_option("-l", "--len", dest="aln_len", type=int, default=50000, help="alignment length (default 50000)"),\
			make_option("-i", "--index", dest="di_index", default="Entropy", help="use [Simpson|Entropy] index, default Entropy"),\
			make_option("-d", "--di-file", dest="di_filename", default=None, help="write out DI to file"),\
			make_option("-o", "--output", dest="tree_filename", default=sys.stderr, help="write out tree to file")
			])

	options, args = parser.parse_args(args=args)
	if options.ecoli_range is not None:
		ecoli_range = eval(options.ecoli_range)  # note: 1-based
	
	if options.di_index not in ("Entropy", "Simpson"):
		print >> sys.stderr, "-i (--index) must either be 'Entropy' or 'Simpson'. Abort!"
		sys.exit(-1)
	
	print >> sys.stderr, "alignment length is", options.aln_len

	df_list = [df for df in DFReader(open(options.df_filename))]
	print >> sys.stderr, "finished reading DF file", options.df_filename

	mask = np.zeros(options.aln_len, dtype=np.float)

	if options.ecoli_only:
		mask[SILVA.Ecoli1542_SILVA100] = 1. # this sets to using ONLY E.coli positions
	else:
		print >> sys.stderr, "Using all positions. Ignore ecoli range."
		mask[:] = 1.

	if options.ecoli_only: 
		mask_lo = ecoli_map[ecoli_range[0]-1]
		mask_hi = ecoli_map[ecoli_range[1]-1]
		mask[:mask_lo] = 0.
		mask[mask_hi+1:] = 0.
		print >> sys.stderr, "taking only positions from E.coli {0}({1})-{2}({3})".format(\
				ecoli_range[0], mask_lo, ecoli_range[1], mask_hi)
#	else:
#		from Solexa_settings import L2
#		ecoli_map = filter(lambda i: L2[i]%1>=0 and 358 <= L2[i] <= 514, xrange(520))
#		import temp_utils
#		nzs = temp_utils.create_threshold_mask_for_df_list(df_list, 1000)
#		raw_input(nzs)
#		mask[:] = 0
#		mask[nzs] = 1.

	V_ecoli = mask.nonzero()[0] # this is the mask we will be using

	for df in df_list:
		#print >> sys.stderr, "changing vec mask for", df.name
		df.change_vec_mask(V_ecoli)
	c = Cluster(df_list, method=options.di_index, threshold=0)
	if options.di_filename is not None:
		print >> sys.stderr, "writing DI to", options.di_filename
		c.write_DI(options.di_filename)
	c.run_till_end()

	t = c.trees[0]
	with open(options.tree_filename, 'w') as f:
		f.write(str(t) + '\n')
	print >> sys.stderr, "output tree written to", options.tree_filename
	return t

if __name__ == "__main__":
	main()

