import os
import sys
import glob
import random
from collections import defaultdict

DIR = '/home/etseng/silo/FH_Meredith/CODE/DigFinger/'
R_SCRIPT_HCLUST = os.path.join(DIR, '20100727_hclust_then_write.R')

def read_sample_rarefaction_file(filename):
	"""
	The file should be:
	<sample-size/real>, <iteration>, <comma-separated list of DI>
	Read the file and return as a dictionary of size --> list of DI strings
	"""
	d = defaultdict(lambda: [])
	with open(filename) as f:
		for line in f:
			size, iter, di = line.strip().split(',', 2)
			if iter == 'real':
				d['real'].append( di )
			else:
				d[int(size)].append( di )
	return d

def jackknifing(d_per_sample, samples, size, output_filename):
	"""
	Randomly pick a DI from size for each of the sample
	"""
	f = open( os.tempnam(), 'w' )
	for sample in samples: 
		print >> sys.stderr, "sample is", sample
		mask = ",".join(map(str,random.choice( d_per_sample[sample][size] ).split(',')[357:514])) # TODO: delete later, testing just using 28-428 (for 27F pyro)
		f.write("{0},{1}\n".format(sample, mask))
	f.close()
	os.system("R --save --input={0} --output={1} < {2}".format( f.name, output_filename, R_SCRIPT_HCLUST ))

def read_hclust_merge_file(filename, num_lines):
	jack_count = defaultdict(lambda: 0) # order tuple --> count
	with open(filename) as f:
		while 1:
			line = f.readline().strip() # this is the header, ignore it unless it's the end
			if len(line) == 0: break
			nodes = {}
			for i in xrange(num_lines):
				iter, a, b = map(int, f.readline().strip().replace('\"', '').split()) # a and b were merged at this point
				print iter,a,b
				if a > 0: # non-singleton merging
					a = nodes[a]
				else: # singleton merging
					a = [-a]
				if b > 0: b = nodes[b]
				else: b = [-b]
				nodes[iter] = a + b
			# now that we've finished reading one permutation's hclust result
			# count how many times the same node appeared
			# remember that the key has to be an ordered TUPLE
			for lst in nodes.itervalues():
				lst.sort()
				jack_count[ tuple(lst) ] += 1
	return dict(jack_count)

if __name__ == "__main__":
	num_permutations = 100
	sample_sizes = [100,500,1000,2000,2500,3000,3500]#[10,50,100,500,1000,2000,3000,4000,5000,6000,7000,8000,9000]

	pattern = sys.argv[1] #pattern = '*4*-*.Entropy.sampled.rarefaction'
	output_filename = sys.argv[2]

	d_per_sample = {}
	for file in glob.iglob(pattern):
		print >> sys.stderr, "reading sample file {0}".format(file)
		d_per_sample[file] = read_sample_rarefaction_file(file)

	sample_names = d_per_sample.keys()
	sample_names.sort() # we keep this sorted list so the cluster labels are always in the same order
	num_samples = len(d_per_sample)

	tmp_filename = os.tempnam()
	jackknifing(d_per_sample, sample_names, 'real', tmp_filename)
	real_nodes = read_hclust_merge_file(tmp_filename, num_samples-1).keys()
	real_nodes.sort(key=lambda x: len(x))
	os.system("rm " + tmp_filename)
	del tmp_filename

	result = defaultdict(lambda: []) # node --> jack_count (in percentage) ex: (1,2) --> .99
	for sample_size in sample_sizes:
		tmp_filename = os.tempnam()
		for p in xrange(num_permutations):
			jackknifing(d_per_sample, sample_names, sample_size, tmp_filename)
		jack_count = read_hclust_merge_file(tmp_filename, num_samples-1)
		os.system("rm " + tmp_filename)
		# node i always refers to sample_names[i-1]
		early_stop = True
		for n in real_nodes:
			try:
				result[n].append(jack_count[n]*1./num_permutations )
			except KeyError: # this means there's not even a single instance of this branch!
				result[n].append(0.)

	with open(output_filename, 'w') as f:
		f.write("NODE\t" + "\t".join(map(str, sample_sizes)) + '\n')
		for n in real_nodes:
			# node i refers to sample_names[i-1]
			node_name = '(' + ",".join(sample_names[i-1] for i in n) + ')'
			f.write(node_name + '\t' + "\t".join(map(str, result[n])) + '\n')

