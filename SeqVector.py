import itertools
#from scipy.sparse import lil_matrix
import numpy as np
import miscIUPAC

class SeqVector:
	mapping = {'N':0, 'A':1, 'T':2, 'G':3, 'C':4}
	rev_mapping = ['N', 'A', 'T', 'G', 'C']
	gap_markers = ['.', '-'] # could be overwritten in settings.py

	def __init__(self, start, end_1, position):
		"""
		start --- now defunct, just set to 0
		end_1 --- this is now the same as len
		position --- should be the position of the first read, as this is
		(defunt)     used later to quickly identify potential mergings
		len --- the length of the reassembled sequence 
		vec --- a 5 x len numpy array where the entry (i, j) denotes the
		        frequency of nucleotide rev_mapping[i] at position j
		conss --- the consensus sequence, is not actually created until the
		          first time consensus_seq() is called
		conss_dirty_bit --- dirty bit flag for consensus_seq()
		read_count --- number of reads used to build this SeqVector, 
		               is automatically incremented every time 
					   add_read or add_seq is called
		"""
		self.start = start
		self.end_1   = end_1  # +1 after the real end
		self.position= position
		self.len   = end_1 - start
		self.vec   = np.zeros((5, self.len), dtype=np.int)
		self.conss = ''
		self.conss_dirty_bit = True
		self.read_count = 0
		self.ids_at_position = []

	def change_vec_mask(self, mask):
		self.end_1 = self.len = len(mask)
		self.vec = self.vec[:, mask]
	
	def add_read(self, read_obj, start_at_offset, copy=1):
		"""
		Add a Read obj at <offset> with abundance <copy>
		"""
		self.ids_at_position.append( (start_at_offset, read_obj.id) )
		self.read_count += copy
		rna_seq = read_obj.read.data.replace('T','U')
		for offset, x in enumerate( rna_seq ):
			if x in SeqVector.gap_markers:
				continue
			nt_s = miscIUPAC.IUPAC_DNA_MAP[ x ]
			w = copy*1. / len(nt_s)
			for nt in nt_s:
				try:
					self.vec[SeqVector.mapping[nt], start_at_offset+offset] += w
				except IndexError:
					print >> sys.stderr, "{0} with start_at_offset {1} offset {2} \
							exceeds length of {3}".format(read_obj.id, start_at_offset, offset, self.end_1)
					print >> sys.stderr, read_obj.gapped_read

		self.conss_dirty_bit = True # remember to set the dirty bit

	def add_seq(self, seq, start_at_offset, copy=1):
		"""
		Same effect as add_read, except that it's a seq
		"""
		self.read_count += copy
		for i,x in enumerate(seq):
			if x in SeqVector.gap_markers:
				continue
			nt_s = miscIUPAC.IUPAC_DNA_MAP[ x ]
			w = copy*1. / len(nt_s)
			for nt in nt_s:
				try:
					self.vec[SeqVector.mapping[nt], start_at_offset+i] += w
				except IndexError:
					pass
		self.conss_dirty_bit = True

	def add_to_vec(self, nt, positions, counts):
		"""
		Add to nt counts[i] for every positions[i]
		"""
		if nt not in SeqVector.gap_markers: 
			self.vec[SeqVector.mapping[nt], positions] += counts # this should be faster than the 2-line version below
#			for i,pos in enumerate(positions):
#				self.vec[SeqVector.mapping[nt], pos] += counts[i]

	def remove_seq(self, seq, start_at_offset, copy=1):
		self.read_count -= copy
		for i,x in enumerate(seq):
			if x in SeqVector.gap_markers:
				continue
			nt_s = miscIUPAC.IUPAC_DNA_MAP[ x ]
			w = copy*1. / len(nt_s)
			for nt in nt_s:
				try:
					self.vec[SeqVector.mapping[nt], start_at_offset+i] -= w
				except IndexError:
					pass
		self.conss_dirty_bit = True

	def get_compressed_vec(self, ignoreN):
		"""
		Returns a 1D (length=self.len) vec that collapses 
		all nucleotide counts per position
		Useful for looking at per-position depth of coverage.
		"""
		result = np.zeros(self.len, dtype=np.int)
		for i in xrange(self.len):
			result[i] = sum(self.vec[:, i])
			if ignoreN: result[i] -= self.vec[SeqVector.mapping['N'], i]
		return result

	def consensus_seq(self, gapped=False, percent_cutoff=.7):
		"""
		At each column, use code that represents >=70% of the counts
		"""
		if not self.conss_dirty_bit:
			return self.conss

		seq = ''
		for col in xrange(self.len):
			total_count = self.vec[:, col].sum()
			nz = self.vec[:, col].nonzero()[0]
			if len(nz) == 0 and gapped:
				seq += '-'
			for size in xrange(1,len(nz)+1):
				done = False
				for st in itertools.combinations(nz, size):
					if sum( self.vec[x, col] for x in st ) >= total_count*percent_cutoff:
						seq += miscIUPAC.get_IUPAC_DNA_code( [SeqVector.rev_mapping[i] for i in st] )
						done = True
						break
				if done:
					break
		self.conss = seq
		self.conss_dirty_bit = False
		return seq

