import os
import sys
import glob
import SILVA
from DiversityIndex import DiversityIndexRunner
from Pyro import Pyro, PyroSampler
from utils.versatile_open import versatile_open_func_wrap

ecoli_map = SILVA.Ecoli1542_SILVA100
subsample_sizes = [10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120]
di_runner = DiversityIndexRunner()

@versatile_open_func_wrap([0],['filename'])
def main(filename, options, f):
	print >> sys.stderr, "processing", filename
	pyro = Pyro(os.path.basename(filename), filename, rename_dups=options.rename_dups)
	sampler = PyroSampler(pyro)
	for se in subsample_sizes:
		if se > pyro.nseq:
			break
		for iter in xrange(options.iterations):
			print >> sys.stderr, "subsample size {0}, iter {1}".format(se, iter)
			df = sampler.subsample(se)
			di = di_runner.run(df, method='Simpson', threshold=0)
			f.write("{name},{size},{di}\n".format(\
					name=pyro.name,\
					size=se,\
					di=",".join(str(di[i]) for i in ecoli_map[ecoli_lo:ecoli_hi])))
	# calc & write the real di
	df = pyro.make_DF()
	di = di_runner.run(df, method='Simpson', threshold=0)
	f.write("{name},real,{di}\n".format(\
			name=pyro.name,\
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

	f = open(output, 'w')
	for filename in glob.iglob(pattern):
		main(filename, options, f)
	f.close()

