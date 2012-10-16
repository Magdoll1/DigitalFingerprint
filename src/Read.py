import os, sys
import cPickle
from Bio import SeqIO
import SILVA
EcoliMap = SILVA.Ecoli1542_SILVA100
from utils import versatile_open
open = versatile_open.versatile_open

def remove_low_quality_for_matched(matches, read_count, phreds, min_phred_score, ditched_f=None):
	"""
	Remove any matches (and it's entries from {matches}, {read_count} and {phreds})
	that have 1 or more positions with a superPhred score < {min_phred_score}

	Returns 
	  count         -- total number of reads removed
	  count_unique  -- total number of unique reads removed
	"""
	count = count_unique = 0
	kk = matches.keys()
	for k in kk:
		m = matches[k]
		if any( x < min_phred_score for x in phreds[m.read.tostring()] ):
			count += read_count[m.read.tostring()]
			count_unique += 1
			if ditched_f is not None:
				ditched_f.write("@{id}\n{seq}\n+{id}\n{qual}\n".format( id=k, seq=m.read, \
					qual=m.quality ))
			del matches[k]
			del read_count[m.read.tostring()]
			del phreds[m.read.tostring()]
	return count, count_unique

class NotEcoliPositionError(Exception):
	def __init__(self, value):
		self.value = value

class RefMap:
	def __init__(self, gap_map_filename, aln_length, fasta_filename=None):
		"""
		gap_map_filename --- which be of format per line:
							{ref_seq_id}\t{comma-separated list of ungapped-to-gapped positions}

		aln_length --- original gapped alignment length, must be given
		"""
		self.gap_map_filename = gap_map_filename
		self.gap_map = {}
		self.aln_length = aln_length # full gapped alignment length, should be CONSTANT for all ref seqs!
		self.read_gap_map()
		if fasta_filename is not None:
			self.fasta = SeqIO.to_dict(SeqIO.parse(open(fasta_filename), 'fasta'))
			self.fasta_filename = fasta_filename

	def read_gap_map(self):
		f = open(self.gap_map_filename)
		for line in f:
			ref_seq_id, positions = line.strip().split('\t')
			positions = map(int, positions.split(','))
			assert all(map(lambda p: 0 <= p < self.aln_length, positions))
			self.gap_map[ref_seq_id] = positions
		f.close()

	def ungapped_to_gapped(self, ref_seq_id, pos):
		return self.gap_map[ref_seq_id][pos]

	def ungapped_to_ecoli(self, ref_seq_id, pos):
		gapped_pos = self.ungapped_to_gapped(ref_seq_id, pos)
		try:
			return EcoliMap.index(gapped_pos)
		except ValueError:
			raise NotEcoliPositionError, "Not an Ecoli position, this needs to be CATCHED!"

class Read:
	def __init__(self, id, seq, copy=1, ref_seq_id=None, offset=None, phred=None):
		self.id = id
		self.seq = seq
		self.phred = phred
		self.ref_seq_id = ref_seq_id
		self.offset = int(offset)
		self.copy = copy

from DF import DF
class ReadDF(DF):
	"""
	Inherits DF to handle the extra (BowTie) read which needs RefMap to add to DF properly
	"""
	def __init__(self, name, refmap, *args):
		DF.__init__(self, name, refmap.aln_length, *args)
		self.refmap = refmap

	def add_read_to_vec(self, read, copy=None):
		"""
		read is a Read object, if copy is None, then read.copy is used
		"""
		for i,s in enumerate(read.seq):
			# the i-th non-gapped position for ref_seq_id starting at offset read.offset
			gapped_pos = self.refmap.ungapped_to_gapped(read.ref_seq_id, read.offset + i)
			DF.add_to_vec(self, nt=s, positions=[gapped_pos], counts=[read.copy if copy is None else copy])

	def add_read_to_vec_using_ref(self, read):
		"""
		match is a BowTieMatch
		instead of adding the match's seq itself...use the ref seq >____<
		"""
		i = read.offset
		for p in self.refmap.gap_map[read.ref_seq_id][read.offset:(read.offset+len(read.seq))]:
			s = self.refmap.fasta[read.ref_seq_id].seq[i]
			if s=='U': s='T'
			if s not in ('A','T','C','G'): s='N'
			DF.add_to_vec(self, nt=s, positions=[p], counts=[read.copy])
			i += 1
			
from Bio.Seq import Seq
from Solexa_settings import MIN_PHRED_SCORE, BOWTIE_PHRED_OFFSET, BowTieMatch
import numpy as np
class ReadsDict:
	def __init__(self, refmap):
		self.M = {} # position --> list of Read objects
		self.refmap = refmap
		self.filenames = []

	def __iter__(self):
		for read_list in self.M.itervalues():
			for read in read_list:
				yield read

	@property
	def count(self):
		"""
		Counts the total number of reads in M
		"""
		return sum(read.copy for read in self.__iter__())

	def special_read(self, bowtie_filename, inhouse_pattern):
		import Solexa_process as sp
		print >> sys.stderr, "readin bowtie"
		matches2, read_count2, phreds2 = sp.gather_reads_BowTie(bowtie_filename)
		sp.remove_low_quality_for_matched(matches2, read_count2, phreds2, MIN_PHRED_SCORE, sys.stderr)
		print >> sys.stderr, "readin inhouse"
		matches, read_count, phreds, unaligned = sp.Shelp.compress_identical_reads(inhouse_pattern)
		sp.remove_low_quality_for_matched(matches, read_count, phreds, MIN_PHRED_SCORE, sys.stderr)
		

	def read_bowtie_match_pickle(self, filename):
		"""
		<filename> should be a pickle that unpickles to:
 		   matches, read_count, phreds
		matches should be a list of (<seq id>, <bowtie match>)
		and read_count should be a dict of seq --> count
		"""
		self.filenames.append(filename)
		with open(filename) as f:
			matches, read_count, phreds = cPickle.load(f)
		print >> sys.stderr, "removing low quality reads with score < {0}".format(MIN_PHRED_SCORE)
		remove_low_quality_for_matched(matches, read_count, phreds, MIN_PHRED_SCORE, None)
		for id, m in matches.iteritems():
			gapped_pos = self.refmap.ungapped_to_gapped(m.ref_seq_id, m.offset)
			if gapped_pos not in self.M:
				self.M[gapped_pos] = []
			read = Read(id, seq=m.read.tostring(), ref_seq_id=m.ref_seq_id, offset=m.offset, \
					copy=read_count[m.read.tostring()])
			self.M[gapped_pos].append(read)

	def read_bowtie_output(self, filename):
		"""
		Modeled after my old Solexa_process.gather_reads_BowTie

		BowTie output file having fields:
		0. read name     1. strand     2. ref name
		3. 0-based offset (w.r.t + strand)     4.seq     5.quality string
		7. # of reserved (ignore in most cases)   7. mismatch descriptors

	    See http://bowtie-bio.sourceforge.net/manual.shtml#algn_out for output format.
		"""
		self.filenames.append(filename)
		matches, read_count, phreds = {}, {}, {}
		f = open(filename)
		for line in f:
			raw = line.strip().split()
			if len(raw) == 5:
				id, strand, ref_seq_id, offset, seq = raw
				qual = [BOWTIE_PHRED_OFFSET] * len(seq) # pretend perfect quality
			else:
				id, strand, ref_seq_id, offset, seq, qual = raw[:6]
				qual = [ord(x) - BOWTIE_PHRED_OFFSET for x in qual]
			if seq in read_count:
				read_count[seq] += 1
				phreds[seq] += qual
			else:
				read_count[seq] = 1
				phreds[seq] = np.array(qual)
				matches[id] = BowTieMatch(id, strand, ref_seq_id, int(offset), Seq(seq), None, None)
		print >> sys.stderr, "removing low quality reads with score < {0}".format(MIN_PHRED_SCORE)
		remove_low_quality_for_matched(matches, read_count, phreds, MIN_PHRED_SCORE, None)
		for id, m in matches.iteritems():
			gapped_pos = self.refmap.ungapped_to_gapped(m.ref_seq_id, m.offset)
			if gapped_pos not in self.M:
				self.M[gapped_pos] = []
			read = Read(id, seq=m.read.tostring(), ref_seq_id=m.ref_seq_id, offset=m.offset, \
					copy=read_count[m.read.tostring()])
			self.M[gapped_pos].append(read)

import bisect
import random
import functools
from Sampler import Sampler
class ReadsSampler(Sampler):
	def __init__(self, readdict):
		self.readdict = readdict
		# whole_pos_index is a sorted list of positions from M
		self.whole_pos_index = self.readdict.M.keys()
		self.whole_pos_index.sort()
		# whole_pos_acc[i] is the cumulative # of read(with identical copies considered) up to i-th position
		# as indicated by whole_pos_index[i]
		M = self.readdict.M
		k0 = self.whole_pos_index[0]
		self.whole_pos_acc = [sum(read.copy for read in M[k0])]
		# per_pos_acc is a dict of position --> ordered cumulative # of copies per M[position]
		self.per_pos_acc = {k0: [M[k0][0].copy]}
		for read in M[k0][1:]:
			self.per_pos_acc[k0].append(self.per_pos_acc[k0][-1] + read.copy)
		# now fill up the rest 
		for k in self.whole_pos_index[1:]:
			self.per_pos_acc[k] = [M[k][0].copy]
			for read in M[k][1:]:
				self.per_pos_acc[k].append(self.per_pos_acc[k][-1] + read.copy)
			self.whole_pos_acc.append(self.whole_pos_acc[-1] + self.per_pos_acc[k][-1] )

		self.n = self.whole_pos_acc[-1] # total number of reads
		self.partial_sampling_func = functools.partial(ReadsSampler.sampling_func,\
				M, self.whole_pos_index, self.whole_pos_acc, self.per_pos_acc)

	@staticmethod
	def sampling_func(M, whole_pos_index, whole_pos_acc, per_pos_acc, i):
		"""
		Finds the i-th read (identical copies considered) in M
		Returns (position of i-th read, i-th read)
		"""
		ind = bisect.bisect_left(whole_pos_acc, i)
		pos = whole_pos_index[ind]
		if ind > 0: i -= whole_pos_acc[ind-1]
		return pos, M[pos][bisect.bisect_left(per_pos_acc[pos], i)]

	def subsample(self, se):
		"""
		Requires 
		partial func --- gives the i-th item (pos, Read obj) in M
		n --- total number of seqs in C
		se --- # of samples to draw
		"""
		df = ReadDF('noname', self.readdict.refmap)
		for i in random.sample(xrange(1, self.n+1), min(se, self.n)):
			pos, read = self.partial_sampling_func(i)
			df.add_read_to_vec(read,copy=1) # important to remember to use just this ONE copy!!!
		return df
