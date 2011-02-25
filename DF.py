import os
import sys
import math
import numpy as np
from SeqVector import SeqVector

class DF(SeqVector):
	number_of_nucleotides = len(SeqVector.mapping)

	def __init__(self, name, len, vec=None):
		self.name = name # unique identifier
		self.len = len # length of vector
		self.__seqvector = SeqVector(0, self.len, 0)
		if vec is not None:
			self.__seqvector.vec = vec
		self.annotations = {}

	def __add__(self, other):
		"""
		Adding one DF to self
		"""
		self.__seqvector.vec += other.__seqvector.vec
		return self

	def __iter__(self):
		"""
		Iterate through the nucleotides
		"""
		for nt in SeqVector.rev_mapping:
			yield nt

	def __len__(self):
		return self.len

	def __str__(self):
		return self.name + '\n' + \
				str(self.__seqvector.vec)

	def __getitem__(self, nt):
		return self.__seqvector.vec[SeqVector.mapping[nt], :]

	def add_annotation(self, k, v):
		self.annotations[k] = v

	def add_to_vec(self, nt, positions, counts):
		"""
		*add* to vec[nt] every counts[i] for positions[i]
		"""
		self.__seqvector.add_to_vec(nt, positions, counts)

	def change_vec_mask(self, mask):
		self.len = len(mask)
		self.__seqvector.change_vec_mask(mask)

	def consensus_seq(self, gapped=False, percent_cutoff=.7):
		return self.__seqvector.consensus_seq(gapped, percent_cutoff)

	def get_compressed_vec(self, ignoreN=True):
		"""
		Calls self.__seqvector.get_compressed_vec
		"""
		return self.__seqvector.get_compressed_vec(ignoreN)

	def get_counts_at_pos(self, i, ignoreN=True):
		"""
		Return count of position i (slice)
		"""
		result = self.__seqvector.vec[:, i]
		if ignoreN:
			nt_ind = SeqVector.mapping['N']
			inds = [i for i in xrange(DF.number_of_nucleotides) if i!=nt_ind]
			return result[inds]
		else:
			return result

	def get_vec_diff_sqsum(self, other, ignoreN=True):
		"""
		Returns the sum of {squared of sum diffs} of two DFs
		"""
		#t = self.vec*1./self.get_compressed_vec(ignoreN) - other.vec*1./other.get_compressed_vec(ignoreN)
		t = self.__seqvector.vec - other.__seqvector.vec
		nt_ind = SeqVector.mapping['N'] if ignoreN else -1
		result = 0
		for j in xrange(self.len):
			c = math.sqrt(sum(t[i,j]**2 for i in xrange(DF.number_of_nucleotides) if i!=nt_ind))
			if not np.isnan(c): result += c
		return result

	def normalized_vec(self, ignoreN=True):
		p = self.get_compressed_vec(ignoreN)
		for i in xrange(self.len):
			p[i] = max(1, p[i])
		self.__seqvector.vec = self.__seqvector.vec * 1. / p
		# the code below was my feeble attempt to correct for the fact
		# that now each positions counts may not add up to 1....*sigh*
#		ind = SeqVector.mapping['A']
#		n_ind = SeqVector.mapping['N']
#		for i in xrange(self.len):
#			if ignoreN:
#				p = 1. - self.__seqvector.vec[:, i].sum() + self.__seqvector.vec[n_ind, i]
#			else:
#				p = 1. - self.__seqvector.vec[:, i].sum()
#			self.__seqvector.vec[ind, i] += p

	def normalized_vec_add(self, other, vec_pre_normalized, ignoreN):
		if vec_pre_normalized:
			self.__seqvector.vec = self.__seqvector.vec * 0.5 + \
					other.__seqvector.vec * 0.5
		else:
			self.__seqvector.vec = self.__seqvector.vec * 0.5 / self.get_compressed_vec(ignoreN)
			self.__seqvector.vec += other.__seqvector.vec * 0.5 / other.get_compressed_vec(ignoreN)

	@property
	def nonzero(self):
		"""
		Return a list of non-zero columns
		because .vec is row(nucleotide) x column(position)
		we take the unified set of columns and return them
		"""
		_x = self.__seqvector.vec.nonzero()[1]
		_x = list(set(_x)) # uniquify them
		_x.sort() # sort positions
		return _x

	@staticmethod
	def nucleotides():
		for nt in SeqVector.mapping: yield nt

	@property
	def vec(self):
		return self.__seqvector.vec

	@property
	def nt_count(self, ignoreN=True):
		inds = range(self.number_of_nucleotides)
		if ignoreN:
			inds.remove(SeqVector.mapping['N'])
		return self.__seqvector.vec[inds, :].sum()


	def assign_vec(self, vec):
		self.__seqvector.vec = vec


class DFReader:
	def __init__(self, f, delimiter=','):
		self.f = f
		self.delimiter = delimiter

	def __iter__(self):
		while 1:
			yield self.next()

	def next(self):
		name = self.f.readline().strip()
		if len(name) == 0:
			raise StopIteration, "read a blank in the file. EOF!"
		length = int(self.f.readline().strip())

		# remember to make the positions 0-based from 1-based
		positions = map(lambda i: int(i)-1, self.f.readline().strip().split(self.delimiter))
		len_positions = len(positions)

		# make a new DF object
		df = DF(name,length)

		for i in xrange(len(SeqVector.rev_mapping)):
			raw = self.f.readline().strip().split(self.delimiter)
			if len(raw) != len_positions + 1:
				raise ValueError, "Each nucleotide line for sample {0} should have exactly {1}+1 elements!".format(\
						name, len_positions)
			nucleotide = raw[0]
			counts_vec = map(int, raw[1:])
			# must add each of counts_vec[i] to positions[i] for the designated nucleotide
			df.add_to_vec(nucleotide, positions, counts_vec)

		# read the annotations if there are any
		while True:
			cur = self.f.tell()
			line = self.f.readline().strip()
			if len(line) == 0:
				break
			if line.startswith('#ANN='):
				k = line[len('#ANN='):]
				line = self.f.readline().strip()
				assert(line.startswith('#ANN:'))
				v = line[len('#ANN:'):]
				df.add_annotation(k, v)
			else:
				self.f.seek(cur)
				break
		return df

class DFWriter:
	def __init__(self, f, delimiter=',', concise=True):
		self.f = f
		self.delimiter = delimiter
		self.concise = concise # concise means only non-zero positions will be written

	def write(self, df):
		assert isinstance(df, DF)

		self.f.write(df.name)
		self.f.write('\n')
		self.f.write(str(df.len))
		self.f.write('\n')

		# remember positions are actually 0-based, but when written out should be 1-based
		if self.concise:
			positions = df.nonzero
		else:
			positions = range(df.len) 
		self.f.write(self.delimiter.join(map(lambda i: str(i+1), positions))) # write as 1-based
		self.f.write('\n')

		for nucleotide in SeqVector.rev_mapping:
			self.f.write(nucleotide + self.delimiter)
			self.f.write(self.delimiter.join(str(df[nucleotide][i]) for i in positions))
			self.f.write('\n')

		for k, v in df.annotations.iteritems():
			self.f.write("#ANN={0}\n".format(k))
			self.f.write("#ANN:{0}\n".format(v))
	
	def writes(self, df_list):
		for df in df_list:
			self.write(df)

def test():
	r = DFReader(open('test'))
	df = r.next()
	print df
	print r.next()
#test()

