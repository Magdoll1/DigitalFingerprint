import os
import sys
import glob
import SILVA
from DiversityIndex import DiversityIndexRunner
from Read import Read, ReadsDict, ReadsSampler, RefMap, ReadDF
from utils.versatile_open import versatile_open_func_wrap

ecoli_map = SILVA.Ecoli1542_SILVA100
subsample_sizes = [10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120]
di_runner = DiversityIndexRunner()

@versatile_open_func_wrap([0],['filename'])
def main(filename, options, f, refmap, ecoli_lo, ecoli_hi):
	from Read import ReadDF
	print >> sys.stderr, "processing", filename
	readdict = ReadsDict(refmap)
	readdict.read_bowtie_output(filename)
	# need to remove all reads that are NOT within the E.coli range
	to_del = filter(lambda x: x < ecoli_map[ecoli_lo] or x >= ecoli_map[ecoli_hi], readdict.M)
	for x in to_del:
		del readdict.M[x]
	print >> sys.stderr, "now {0} reads left".format(readdict.count)

	sampler = ReadsSampler(readdict)

	for se in subsample_sizes:
		for iter in xrange(options.iterations):
			print >> sys.stderr, "subsample size {0}, iter {1}".format(se, iter)
			df = sampler.subsample(se)
			di = di_runner.run(df, method='Simpson', threshold=0)
			f.write("{name},{size},{di}\n".format(\
					name=filename,\
					size=se,\
					di=",".join(str(di[i]) for i in ecoli_map[ecoli_lo:ecoli_hi])))
	# calc & write the real di
	df = ReadDF(filename, refmap)
	for read in readdict:
		df.add_read_to_vec(read)
	di = di_runner.run(df, method='Simpson', threshold=0)
	f.write("{name},real,{di}\n".format(\
			name=filename,\
			di=",".join(str(di[i]) for i in ecoli_map[ecoli_lo:ecoli_hi])))

if __name__ == "__main__":
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option("-p", "--pattern", dest="pattern", help="pattern for list of input files")
	parser.add_option("-o", "--output", dest="output", help="output filename")
	parser.add_option("--dup-ok", dest="rename_dups", action="store_true", default=False, help="ok to have duplicate IDs (renames them)")
	parser.add_option("-i", "--iterations", dest="iterations", type=int, help="number of iterations to run subsampling on")
	parser.add_option("-r", "--ecoli-range", dest="ecoli_range", help="Ecoli range to run subsampling/DI on")

	options, args = parser.parse_args()
	pattern, output = options.pattern, options.output

	ecoli_lo, ecoli_hi = eval(options.ecoli_range) # NOTE: 1-based, must convert to 0-based for ecoli_map
	ecoli_lo = int(ecoli_lo)- 1
	ecoli_hi = int(ecoli_hi)

	print >> sys.stderr, "reading refmap....."
	refmap = RefMap(os.environ['PCODE'] + '/Silva/SILVA104.fece_augmented.fasta.gap_map.bz2', aln_length=50000)

	f = open(output, 'w')
	for filename in glob.iglob(pattern):
		main(filename, options, f)
	f.close()

