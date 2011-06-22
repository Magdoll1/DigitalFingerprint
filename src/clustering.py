import sys
import math
import numpy as np
from newick import tree
from DiversityIndex import DiversityIndexRunner, SimpsonIndex, EntropyIndex
"""
Modified hierarchical clustering using Digital Fingerprinting
"""
class Cluster:
	def __init__(self, df_list, **kwargs):
		if df_list is None:
			print >> sys.stderr, "called with a None, I hope you know what you're doing!"
			print >> sys.stderr, "calling init_from_di_list later perhaps?"
			self.df_list = None
			return
		self.df_list = df_list
		self.original_names = [df.name for df in self.df_list]
		self.mask = kwargs['mask'] if 'mask' in kwargs else 1.
		self.method = kwargs['method'] if 'method' in kwargs else 'Simpson'
		self.threshold = kwargs['threshold'] if 'threshold' in kwargs else 10

		self.m = len(df_list) # number of samples (rows)
		self.n = df_list[0].len # length of the DF vectors (columns)

		self.runner = DiversityIndexRunner(self.mask)

		self.trees = [tree.Leaf(self.df_list[i].name) for i in xrange(self.m)]

		self.X = np.zeros((self.m, self.n), dtype=np.float)
		for i,df in enumerate(self.df_list):
			#print >> sys.stderr, "normalizing {0}....".format(df.name)
			df.normalized_vec(ignoreN=True)
			di = self.runner.run(df, method=self.method, threshold=self.threshold, \
					vec_pre_normalized=True, ignoreN=True)
			self.X[i, :] = di	

		# calculate the initial distance matrix
		self._dist = np.zeros((self.m, self.m), dtype=np.float)
		for i in xrange(self.m):
			self._dist[i, i] = float("inf")
			for j in xrange(i+1, self.m):
				# method 1: Euclidean distance between DIs
				d = math.sqrt(sum(x**2 for x in self.X[i,:]-self.X[j,:]))
				# method 2: sum of sum of distances squared between DFs
				#d = self.df_list[i].get_vec_diff_sqsum(self.df_list[j])
				self._dist[i, j] = d
				self._dist[j, i] = d
	
	def init_from_di_list(self, di_dict, **kwargs):
		"""
		alternative __init__ taking a dict sample_name --> array of DI as input
		if this init is used, we're expecting to run UPGMA clustering
		"""
		self.original_names = di_dict.keys()
		self.original_names.sort()
		self.mask = kwargs['mask'] if 'mask' in kwargs else 1.
		self.method = kwargs['method'] if 'method' in kwargs else 'Simpson'
		self.threshold = kwargs['threshold'] if 'threshold' in kwargs else 10
		self.m = len(di_dict)
		self.n = len(di_dict.itervalues().next())
		self.trees = [tree.Leaf(x) for x in self.original_names]
		self.X = np.zeros((self.m, self.n), dtype=np.float) 
		# fill up X using di_dict
		for i in xrange(self.m):
			self.X[i] = di_dict[self.original_names[i]]
		self._dist = np.zeros((self.m, self.m), dtype=np.float)
		for i in xrange(self.m):
			self._dist[i, i] = float("inf")
			for j in xrange(i+1, self.m):
				# method 1: Euclidean distance between DIs
				d = math.sqrt(sum(x**2 for x in self.X[i,:]-self.X[j,:]))
				# method 2: sum of sum of distances squared between DFs
				#d = self.df_list[i].get_vec_diff_sqsum(self.df_list[j])
				self._dist[i, j] = d
				self._dist[j, i] = d
	
	def write_DI(self, output_filename, mask=None):
		with open(output_filename, 'w') as f:
			for i, name in enumerate(self.original_names):
				di = self.X[i, ] if mask is None else self.X[i, mask]
				f.write(name + ',')
				di.tofile(f, sep=",")
				f.write('\n')
						
	def run_one_cluster_step(self):
		d = self._dist.argmin()
		i, j = d / self.m, d % self.m
		_min_val = self._dist[i, j]
		if _min_val == float("inf"):
			raise StopIteration, "done!"
		#print >> sys.stderr, "combining {0} and {1}".format(self.trees[i], self.trees[j])
		
		# merge j into i
		size_i = len(self.trees[i].get_leaves())
		size_j = len(self.trees[j].get_leaves())
		t = tree.Tree()
		t.add_edge((self.trees[i], 0, _min_val/2)) # (subtree-i, bootstrap=0, branch length=dist)
		t.add_edge((self.trees[j], 0, _min_val/2)) 
		self.trees[i] = t
		self.trees[j] = None
		# NEW!!! instead of just adding df_list[j] to df_list[i], normalize the counts FIRST!!!
		if self.df_list is not None:
			self.df_list[i].normalized_vec_add(self.df_list[j], vec_pre_normalized=True, ignoreN=True)

#			print "before", self.X[i, ]
			self.X[i] = self.runner.run(self.df_list[i], method=self.method, threshold=self.threshold,\
					vec_pre_normalized=True, ignoreN=True)
#			print("merged {0} and {1}".format(i, j))
#			print "new vec is now", self.X[i, ]
		
#		self._dist[j, :] = float("inf")
#		self._dist[:, j] = float("inf")
		for k in xrange(self.m):
			if k==i or k==j or self.trees[k] is None: continue
			# method 1:
			#d = math.sqrt(sum(x**2 for x in self.X[i,:]-self.X[k,:]))
			# method 2:
			#d = self.df_list[i].get_vec_diff_sqsum(self.df_list[k])
			# method 3: UPGMA
			d = (self._dist[k, i] * size_i + self._dist[k, j] * size_j) / (size_i + size_j)
			# method 4: complete linkage
			#d = max(self._dist[k, i], self._dist[k, j])
			#print >> sys.stderr, "using Euclidean dist: {0}, using vecdiff: {1}".format(\
			#		math.sqrt(sum(x**2 for x in self.X[i,:]-self.X[k,:])), d)
			self._dist[i, k] = d
			self._dist[k, i] = d
		self._dist[j, :] = float("inf")
		self._dist[:, j] = float("inf")

	def run_till_end(self):
		while len(self.trees) > 1:
			try:
				self.run_one_cluster_step()
			except StopIteration:
				break

def main(args=None):
	from DF import DFReader
	import SILVA 
	ecoli_map = SILVA.Ecoli1542_SILVA100
	
	from optparse import make_option, OptionParser
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

