from abc import ABCMeta, abstractmethod
from DF import DF, DFReader
import numpy as np
import math

class DiversityIndex:
	"""
	Abstract class
	"""
	__metaclass__ = ABCMeta

	@abstractmethod
	def calc_index(self, df, threshold):
		raise NotImplemented, "Should not call the abstract method calc_index!"

class SimpsonIndex(DiversityIndex):
	@staticmethod
	def calc_index(df, threshold, vec_pre_normalized, ignoreN):
		D = np.zeros(df.len, dtype=np.double)
		for i in xrange(df.len):
			_x = df.get_counts_at_pos(i, ignoreN)
			if vec_pre_normalized: # vec is already in percentages
				D[i] = 1 - sum(s**2 for s in _x)
			else:
				n = sum(_x)
				if n >= threshold and n > 1: 
					D[i] = 1 - sum(s * (s - 1) for s in _x) * 1. / (n * (n-1))
		return D

class EntropyIndex(DiversityIndex):
	"""
	Similar to diversity_index except that diversity_index uses Simpson Index
	and this uses the entropy information which is 
		sum_{s\in A,T,C,G} P(s)*log( P(s)/Q(s) )
	where Q is the background nucleotide distribution and
	      P is the current column's nucleotide distribution
	Right now we use log 2.
	(N_cutoff is NOT used!)
	"""
	@staticmethod
	def calc_index(df, threshold):
		D = np.zeros(df.len, dtype=np.double)
		log2 = lambda x: math.log(x, 2)
		# first we calculate the background distribution
		totaln = sum(sum(df[nt]) for nt in df) * 1.
		Q = dict((nt, sum(df[nt])/totaln) for nt in df)

		# now calculate per-column entropy
		for i in xrange(df.len):
			n = sum(df.get_counts_at_pos(i, ignoreN=True)) * 1.
			if n >= threshold and n > 1:
				for nt in df: 
					p_i = df[nt][i] / n
					q_i = Q[nt]
					if p_i > 0:
						D[i] += p_i * log2(p_i / q_i)
		return D

class DiversityIndexRunner:
	runners = {}
	def __init__(self, mask=None):
		self.mask = 1. if mask is None else mask

	@staticmethod
	def register_class_name(_string, _class):
		DiversityIndexRunner.runners[_string] = _class

	def run(self, df, **kwargs):
		options = {'method': 'Simpson',\
				'threshold': 0,\
				'vec_pre_normalized': True,\
				'ignoreN': True}
		options.update(kwargs)
		return DiversityIndexRunner.runners[options['method']].calc_index(df, options['threshold'], options['vec_pre_normalized'], options['ignoreN']) * self.mask

DiversityIndexRunner.register_class_name('Simpson', SimpsonIndex)
DiversityIndexRunner.register_class_name('Entropy', EntropyIndex)


# --- test code below --- #
#df = DFReader(open('SILVA100_justHumanCrap.full.DF')).next()
#test = DiversityIndexRunner()
#di = test.run(df, method='Simpson', threshold=0)
