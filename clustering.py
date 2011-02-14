import sys
import math
import numpy as np
import hcluster #TODO: handle it if hcluster isn't available -- homemade pdist?
from newick import tree
from DiversityIndex import DiversityIndexRunner, SimpsonIndex, EntropyIndex
"""
Modified hierarchical clustering using Digital Fingerprinting
"""
#def dist(di_list):
#	d = []
#	n = len(di_list)
#	for i in xrange(n-1):
#		for j in xrange(i+1, n):
#			d.append(math.sqrt((di_list[i] - di_list[j])**2))
	
def find_index_in_condensed_array_help(cur, length, i):
	if cur < length:
		return (i, cur+i+1)
	else:
		return find_index_in_condensed_array_help(cur-length, length-1, i+1)

def find_index_in_condensed_array(index, length):
	return find_index_in_condensed_array_help(index, length-1, 0)


class Cluster:
	def __init__(self, df_list, **kwargs):
		self.df_list = df_list
		self.original_names = [df.name for df in self.df_list]
		self.mask = kwargs['mask'] if 'mask' in kwargs else 1.
		self.method = kwargs['method'] if 'method' in kwargs else 'Simpson'
		self.threshold = kwargs['threshold'] if 'threshold' in kwargs else 10

		self.m = len(df_list) # number of samples (rows)
		self.n = df_list[0].len # length of the DF vectors (columns)

		self.runner = DiversityIndexRunner(self.mask)

		self.trees = [tree.Leaf(self.df_list[i].name) for i in xrange(self.m)]

		self.X = np.zeros((self.m, self.n), dtype=np.double)
		for i,df in enumerate(self.df_list):
			df.normalized_vec()
			di = self.runner.run(df, method=self.method, threshold=self.threshold)
			self.X[i, :] = di	
			print("self.X[{0}] is {1}".format(i, self.X[i,:]))

		# calculate the initial distance matrix
		#self._dist = hcluster.pdist(self.X, 'euclidean')
		self._dist = np.zeros((self.m, self.m), dtype=np.double)
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
		print >> sys.stderr, "combining {0} and {1}".format(self.trees[i], self.trees[j])
		
		#_min_ind = self._dist.argmin()
		#_min_val = self._dist[_min_ind]
		#i, j = find_index_in_condensed_array(_min_ind, self.m) # guaranteed i < j
		# pop out trees[i], trees[j], merge them into a new one
		size_i = len(self.trees[i].get_leaves())
		size_j = len(self.trees[j].get_leaves())
		t = tree.Tree()
		t.add_edge((self.trees[i], 0, _min_val/2)) # (subtree-i, bootstrap=0, branch length=dist)
		t.add_edge((self.trees[j], 0, _min_val/2)) 
		self.trees[i] = t
		self.trees[j] = None#self.trees.pop(j) # only do this if using pdist array
		# NEW!!! instead of just adding df_list[j] to df_list[i], normalize the counts FIRST!!!
#		self.df_list[i] += self.df_list[j] # OLD way: just add
#		self.df_list[i].normalized_vec()
		self.df_list[i].normalized_vec_add(self.df_list[j])
		# TEMP: just turn it into ratios and average over two
		#self.df_list[i].normalized_vec_add(self.df_list[j])
		#self.df_list.pop(j) # only do this if using pdist-array
		#self.m -= 1 # only do this if using pdist-array

		print "before", self.X[i, ]
		self.X[i] = self.runner.run(self.df_list[i], method=self.method, threshold=self.threshold)
		print("merged {0} and {1}".format(i, j))
		print self.df_list[i].vec[:,-3][1:],
		print self.df_list[j].vec[:,-3][1:]
		print "new vec is now", self.X[i, ]
		print "pos 0", self.X[i,0]
#		raw_input("PRESS ANY KEY")
#		X = self.X.tolist()
#		X.pop(j)
#		self.X = np.array(X)
		# re-run distance calculation
		# TODO: make this NOT dependent on hcluster...!
#		self._dist = hcluster.pdist(self.X, 'euclidean')
		
#		self._dist[j, :] = float("inf")
#		self._dist[:, j] = float("inf")
		for k in xrange(self.m):
			if k==i or k==j or self.trees[k] is None: continue
			# method 1:
			#d = math.sqrt(sum(x**2 for x in self.X[i,:]-self.X[k,:]))
			# method 2:
			#d = self.df_list[i].get_vec_diff_sqsum(self.df_list[k])
			#method 3:
			d = (self._dist[k, i] * size_i + self._dist[k, j] * size_j) / (size_i + size_j)
			print >> sys.stderr, "using Euclidean dist: {0}, using vecdiff: {1}".format(\
					math.sqrt(sum(x**2 for x in self.X[i,:]-self.X[k,:])), d)
			self._dist[i, k] = d
			self._dist[k, i] = d
		self._dist[j, :] = float("inf")
		self._dist[:, j] = float("inf")

		print "dist is:",
		print self._dist
	def run_till_end(self):
		while len(self.trees) > 1:
			try:
				self.run_one_cluster_step()
			except StopIteration:
				break

if __name__ == "__main__":
	import SILVA 
	ecoli_map = SILVA.Ecoli1542_SILVA100
	
	from optparse import make_option, OptionParser
	parser = OptionParser(option_list=[ \
			make_option("-f", "--df-filename", dest="df_filename"), \
			make_option("-r", "--ecoli-range", dest="ecoli_range", help="specify the E.coli position range to use ex:(238,338)"), \
			make_option("-a", "--all-positions", dest="ecoli_only", action="store_false", default=True, help="use all 50,000 positions"),\
			make_option("-d", "--di-file", dest="di_filename", default=None, help="write out DI to file")
			])

	options, args = parser.parse_args()
	if options.ecoli_range is not None:
		ecoli_range = eval(options.ecoli_range)  # note: 1-based

#	mask = np.array(SILVA.DI_Simpson_Ecoli1542_SILVA100)
	from DF import DFReader
	df_list = [df for df in DFReader(open(options.df_filename))]
	print >> sys.stderr, "finished reading DF file", options.df_filename

#	import temp_utils
#	mask = temp_utils.create_threshold_mask_for_df_list(df_list, threshold=100)
	mask = np.zeros(50000, dtype=np.float)#mask = np.zeros(50000, dtype=np.float)

	if options.ecoli_only:
		mask[SILVA.Ecoli1542_SILVA100] = 1. # this sets to using ONLY E.coli positions
	else:
		mask[:] = 1.

	if options.ecoli_only:
		mask_lo = ecoli_map[ecoli_range[0]-1]
		mask_hi = ecoli_map[ecoli_range[1]-1]
		mask[:mask_lo] = 0.
		mask[mask_hi+1:] = 0.
		print >> sys.stderr, "taking only positions from E.coli {0}({1})-{2}({3})".format(\
				ecoli_range[0], mask_lo, ecoli_range[1], mask_hi)
	else:
		from Solexa_settings import L2
		ecoli_map = filter(lambda i: L2[i]%1==0 and 358 <= L2[i] <= 514, xrange(520))
		mask[:] = 0.
		mask[ecoli_map] = 1.
#		mask[:] = 1. # TODO: remove later
		#print >> sys.stderr, "taking all positions from {0}-{1}".format(mask_lo, mask_hi)

#	mask[:6427] = 0. # remove all locations < E.coli 358 (which is 6427)
#	mask[11895:] = 0. # remove all locations > E.coli 514 (which is 11894)
#	mask[:27655] = 0. # remove all locations < E.coli 907
#	mask[34343:] = 0. # remove all locations > E.coli 1073
#	mask[6334:] = 0. # remove all locations > E.coli 338 (which is 6333)
#	mask[:5280] = 0. # remove all locations < E.coli 238 (which is 5280)
#	mask[5280:6334] = 1.
#	mask[:3856] = 0. # remove all locations < E.coli 188 (which is 3856)
#	mask[:2084] = 0. # remove all locations < E.coli 138 (which is 2084)
#	mask[:13151] = 0. # remove all locations < E.coli 528 (which is 13151)
#	mask[26145:] = 0. # remove all locations > E.coli 828 (which is 22059)
#	mask[p] = 1.
#	V2_ecoli_last200bp = mask.nonzero()[0] # this is Ecoli 238-338
#	V6_ecoli = mask.nonzero()[0] # this is Ecoli 907-1073
#	V3_ecoli = mask.nonzero()[0] # this is Ecoli 358-514
#	V4_ecoli = mask.nonzero()[0] # this is Ecoli 528-828
	V_ecoli = mask.nonzero()[0] # this is THE Ecoli mask

	for df in df_list:
#		print >> sys.stderr, "changing vec mask for", df.name
		df.change_vec_mask(V_ecoli)

#	mask[27655:34343] = 1.
#	mask = mask * np.array(SILVA.DI_Simpson_humanCrap_SILVA100)
#	c = Cluster(df_list, method='Simpson', threshold=100, mask=mask)
	c = Cluster(df_list, method='Simpson', threshold=0)
	if options.di_filename is not None:
		print >> sys.stderr, "writing DI to", options.di_filename
		c.write_DI(options.di_filename)

#	for i in xrange(c.m):
#		print(c.original_names[i] + ',' + ",".join(map(str, c.X[i,])))
#	print()
#	d = hcluster.squareform(c._dist)
#	for i in xrange(c.m):
#		print(",".join(map(str, d[i,:])))

	c.run_till_end()
	print c.trees[0]
