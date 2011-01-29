import os
import sys
from Pyro import Pyro, PyroSampler
from clustering import Cluster
class Jackknifer:
	def __init__(self, samplers, mask, method, threshold):
		self.samplers = samplers
		self.mask = mask
		self.method = method
		self.threshold = threshold
		self.trees = {} # subsample size --> list of trees from each run
		
	def run(self, se):
		df_list = [sampler.subsample(se) for sampler in self.samplers]
		c = Cluster(df_list, mask=self.mask, method=self.method, threshold=self.threshold)
		c.run_till_end() # c.trees[0] will be the only tree and the final tree
		if se not in self.trees:
			self.trees[se] = []
		self.trees[se].append(c.trees[0]) 

	def runs(self, se, iter):
		for i in xrange(iter):
			print >> sys.stderr, "iteration {0} for subsampling at size {1}...".format(i, se)
			self.run(se)

if __name__ == "__main__":
	import glob
	import numpy as np
	import SILVA

	mask = np.array(SILVA.DI_Simpson_Ecoli1542_SILVA100)

	samplers = []
	for file in glob.iglob('*.align.quality_filtered'):
		print >> sys.stderr, "reading file {0}....".format(file)
		pyro = Pyro(os.path.basename(file), file)
		pyro_sampler = PyroSampler(pyro)
		samplers.append(pyro_sampler)

	jack = Jackknifer(samplers, mask, 'Simpson', threshold=1000)

	subsample_size = int(sys.argv[1])
	iterations = int(sys.argv[2])
	output_filename = sys.argv[3]
	
	jack.runs(subsample_size, iterations)
	with open(output_filename, 'w') as f:
		for t in jack.trees[subsample_size]:
			f.write(str(t) + '\n')

