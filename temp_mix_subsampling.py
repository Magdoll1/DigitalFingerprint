import os
import sys
import random
from utils.versatile_open import versatile_open_func_wrap
"""
Given a bunch of subsampled files each of which should be of per-line format:
	<name>,<subsample size or 'real'>,<comma-separated list of DI>

Read them, and generate <X> files each of which is one random selection of DI
from each of the subsampled sizes
"""
subsample_sizes = [10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120]

@versatile_open_func_wrap([0],['filename'])
def populate_samples(filename, S, S_names):
	"""
	S is subsample size --> sample name --> list of DIs
	"""
	with open(filename) as f:
		for line in f:
			name,size,di = line.strip().split(',', 2)
			size = int(size) if size!='real' else 'real'
			if size not in S: S[size] = {}
			if name not in S[size]: S[size][name] = []
			S[size][name].append(di)
			S_names.add(name)

def mix_samples(S, S_names, subsample_sizes, output_filename):
	with open(output_filename, 'w') as f:
		for se in subsample_sizes:
			for name in S_names:
				if name in S[se]: 
					x = S[se][name]
				else: # probably becuz this exceeded the largest size, just use 'real'
					x = S['real'][name]
				f.write("{name},{size},{di}\n".format(\
						name=name,\
						size=se,\
						di=random.choice(x)))
		# now write out 'real'
		for name in S_names:
			f.write("{name},real,{di}\n".format(\
					name=name,\
					di=S['real'][name][0]))

if __name__ == "__main__":
	from optparse import OptionParser, make_option

	parser = OptionParser(option_list=[\
			make_option("-f", "--filenames", dest="filenames", help="comma-separated list of input files"),\
			make_option("-i", "--iterations", dest="iterations", type=int),\
			make_option("-o", "--output", dest="output_prefix", help="output prefix")])

	options, args = parser.parse_args()

	S = {}
	S_names = set()
	for filename in options.filenames.split(','):
		populate_samples(filename, S, S_names)

	S_names = list(S_names)
	S_names.sort() # just for astheticism

	for iter in xrange(options.iterations):
		print >> sys.stderr, "mixing iteration {0}....".format(iter)
		mix_samples(S, S_names, subsample_sizes, options.output_prefix+'_iter'+str(iter)+'.txt')

