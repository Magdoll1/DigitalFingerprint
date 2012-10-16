import os,re,sys,bisect
from collections import defaultdict, namedtuple
from Bio import SeqIO
import miscIUPAC

DIR = '/home/etseng/FH_Meredith/SolexaReads/'

ECOLI_HELIX   = DIR + 'refDBs/d.16.b.E.coli.bpseq_1crapEXPANDED'
#CONS_filename = 'obesity_original_CONS_allPhylum.Ecoli_masked.fasta'

filename_refDB_aligned = DIR + 'refDBs/SILVA100_justHumanCrap.1crap_masked.V3region.fna'
#filename_locus_to_org  = DIR + 'refDBs/SILVA100_justHumanCrap.id_to_tax_slv.txt'

GAP_SYMBOLS = ['-', '.', '?']

V3_offset = 0
V3_length = 520 # this is the length of the V3 1crap mask

class HybridCons:
	V3_length = 520
	NO_BP_PRESENT = 1.5
	def __init__(self):
#		self.seq = {}
#		self.pairing = {}
		self.bpseq = read_BPSEQ( ECOLI_HELIX )
		self.mapping = map_position_of_ungapped_to_gapped()
		self.refseqs = SeqIO.to_dict( SeqIO.parse(open(filename_refDB_aligned),'fasta') )
#		self.locus_to_phylum = {}

#		with open(CONS_filename) as handle:
#			for r in SeqIO.parse( handle, 'fasta' ):
#				seq = r.seq.tostring()
#				phylum = r.id[len('CONS_'):]
#
#				self.seq[phylum] = seq
#				print >> sys.stderr, "handling {0}".format(phylum)
#				self.pairing[phylum] = hybrid_structure(self.bpseq, seq)

		# populate locus_to_phylum
#		with open(filename_locus_to_org) as handle:
#			for line in handle:
#				locus,tax = line.strip().split('\t')
#				# phylum is the 2nd level
#				self.locus_to_phylum[locus] = tax.split('/')[1]

	def next_nongap_pos(self, ref_seq_id, i):
		"""
		Find the next non-gap position after (gapped) i-th position of ref_seq_id
		"""
		try:
			ungapped_i = self.mapping[ref_seq_id].index(i)
		except ValueError:
			# find the largest ungapped_i s.t. mapping[ref_seq_id][ungapped_i] < i
			ungapped_i = bisect.bisect(self.mapping[ref_seq_id], i) - 1

		if ungapped_i + 1 < len(self.mapping[ref_seq_id]):
			return self.mapping[ref_seq_id][ungapped_i+1]
		else:
			return None

	def qualify_with_secondary_structure(self, read, ref_seq_id, i, match_size):
		"""
		A read matches to <ref_seq_id> starting at (local V3) position <i>
		Check for each base-pairing within the range of this 51-bp read,
		 whether there's a base-pairing whenever bpseq says there is
		"""
		bped, acc = 0, 0
		mapping = self.mapping[ref_seq_id]
		# the 0-crap (full-length, gapped) position range of this read is thus:
		# mapping[i] + V3_offset -- mapping[i+50] + V3_offset
		try:
			j = i + match_size - 1
			l,h = mapping[i] + V3_offset, mapping[j] + V3_offset
		except IndexError:
			raise Exception, "trying to index {0}[{1}] and {0}[{2}] failed".format(ref_seq_id, i, j)
		for x,y in self.bpseq.iteritems():
			if l <= x <= h and l <= y <= h and x < y:
				try:
					ii = mapping.index( x - V3_offset ) - i
					jj = mapping.index( y - V3_offset ) - i
				except ValueError:
					continue
				# it is possible that ii or jj is beyond the range of read
				# in that case just ignore
				if ii < 0 or jj < 0 or ii >= len(read) or jj >= len(read):
					print >> sys.stderr, "IGNORE some in qualifying 2nd structure cuz out of range"
					continue
				acc += 1
				if miscIUPAC.can_pair(read[ii], read[jj]):
					bped += 1

		if acc == 0: # if there is no bp present, then the ratio counts as 1.5!
			return HybridCons.NO_BP_PRESENT
		else:
			return bped*1. / acc

	def qualify_assembled(self, gapped_seq, offset):
		"""
		Given an assembled, gapped seq that starts at (V3)'s offset-th position
		see how many of the base-pairs match.

		Returns -1 if none of the positions are supposed to be base-paired anyways.
		Otherwise returns a fraction of <should-bp-and-did> / <should-bp>
		"""
		bped, acc = 0, 0
		
		i = V3_offset + offset
		j_1 = i + len(gapped_seq)

		for x,y in self.bpseq.iteritems():
			if i <= x < j_1 and i <= y < j_1 and x < y:
				u = gapped_seq[ x - V3_offset ]
				v = gapped_seq[ y - V3_offset ]
				acc += 1
				if u not in GAP_SYMBOLS and v not in GAP_SYMBOLS and miscIUPAC.can_pair(u, v):
					bped += 1

		if acc == 0:
			return 1.
		else:
			return bped*1. / acc

def map_position_of_ungapped_to_gapped():
	"""
	Returns a dictionary where key is ref_seq_id and value is a list 
	list[i]=j means i-th position in ungapped maps to j-th position in gapped
	"""
	mapping = {} # ref_seq_id --> [ gapped_position for 0,1,2-th position in ungapped ]
	with open(filename_refDB_aligned) as handle:
		for r in SeqIO.parse( handle, 'fasta' ):
			mapping[r.id] = []
			seq2 = r.seq.tostring()
			seq1 = seq2
			for g in GAP_SYMBOLS: seq1 = seq1.replace(g, '')

			j = 0
			for i,x in enumerate(seq1):
				while seq2[j] == '-': j += 1
				assert x == seq2[j]
				mapping[r.id].append( j )
				j += 1

	return mapping

def read_BPSEQ(filename):
	"""
	Reads a BQSEQ file format
	Returns a dict where key --> val are base-paired positions (0-based)

	asumes the first 4 lines are headers, and positions are 0-based in format:

	<i-th>	<nucleotide at i-th>	<j-th, base pairs with i-th>
	"""
	bp_map = {}
	with open(filename) as f:
		# ignore headers
		for x in xrange(4): f.readline()
		for line in f:
			x,nt,y = line.strip().split()
			x = int(x)
			y = int(y)
			bp_map[x] = y
			bp_map[y] = x
	return bp_map

def hybrid_structure(bp_map, seq):
	"""
	Given bp_map (output from read_BPSEQ)
	and seq which is a string (can have gaps)
	returns list L where 
		L[i] is None      if i-th position is gap
		L[i] is (i,nt,0)  if i-th position has no pairing
		L[i] is (i,nt,j)  if i-th position pairs with j-th position
	"""
	seq = seq.upper().replace('T', 'U')
	seq_len = len(seq)
	L = [0] * seq_len

	for i,x in enumerate(seq):
		if x in GAP_SYMBOLS:
			L[i] = None
		elif i not in bp_map:
			L[i] = (i, x, 0)
		else:
			y = seq[ bp_map[i] ]
			if y not in GAP_SYMBOLS and miscIUPAC.can_pair(x, y):
				L[i] = (i, x, bp_map[i])
			else:
				L[i] = (i, x, 0)

	return L

def print_hybrid(L):
	for bb in L:
		if bb is None:
			continue
		i, nt, j = bb
		if j == 0:
			print("{0}\t{1}\t0".format(i+1, nt))
		else:
			print("{0}\t{1}\t{2}".format(i+1, nt, j+1))


if __name__ == "__main__":
#	seq = SeqIO.read( open('obesitu_original_CONS_Actinobacteria.Ecoli_masked.fasta'), 'fasta' ).seq.tostring()
	h = HybridCons()
	ref_seq_id='Unc05j70'
	read='GAACGAGACGCCCTTCGGGGTGTAAAGTTTTGTCAGTGGGGACGAACGAAT'
	print h.qualify_with_secondary_structure(read, ref_seq_id, 68, len(read))
