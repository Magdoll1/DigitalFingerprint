import numpy as np
cimport numpy as np
from Solexa_settings import BOWTIE_PHRED_OFFSET, NT_MAPPING, BowTieMatch
import cPickle

def gather_reads_BowTie(object filename, object refmap, np.ndarray[np.int64_t, ndim=2] seqvec, int phred_cutoff, int min_length, int ecoli_pos_lo, int ecoli_pos_hi):
	"""
	filename --- bowtie output, must have <name>, <strand>, <ref seq id>, <offset>, <read>, <quality>, <the_rest_is_ignored>
	refmap   --- Read.RefMap object
	seqvec   --- the mxn nucleotide count matrix where seqvec[i,j] is count of nucleotide i for position j
	phred_cutoff, min_length --- trim reads up to the first bad base (<phred_cutoff) we see and use it only if the trimmed read
	                             if longer than <min_length>
	"""
	cdef int i, j, pos, ind, len_read
	cdef int discard = 0
	cdef int use = 0
	f = open(filename)
	for line in f:
		name, strand, ref_seq_id, offset, read, quality, junk = line.strip().split('\t', 6)
		len_read = len(read)
		offset = int(offset)
		if refmap.gap_map[ref_seq_id][offset] < ecoli_pos_lo or\
				refmap.gap_map[ref_seq_id][offset] > ecoli_pos_hi:
#			raw_input("discarding read becuz mapped pos " + str(refmap.gap_map[ref_seq_id][offset]))
			discard += 1
			continue
		for i in range(len_read+1):
			if i == len_read:
				break
			phred = ord(quality[i]) - BOWTIE_PHRED_OFFSET
			assert 0 <= phred <= 40
			if phred < phred_cutoff:
				break
		if i >= min_length:
			#print("using the first {0} bases: {1}".format(i, read[:i]))
			for j in range(i):
				pos = refmap.gap_map[ref_seq_id][offset+j]
				nt = read[j]
				if nt not in NT_MAPPING:
					nt = 'N'
				ind = NT_MAPPING[nt]
				#print("adding count to index {0}, pos {1}".format(ind, pos))
				seqvec[ind, pos] += 1
			use += 1
		else:
			#print("discarding {0} becuz i is {1}".format(read, i))
			discard += 1
	f.close()
	return use, discard

def subsample_reads_BowTie_prepare(object filename, object refmap, int phred_cutoff, int min_length, int ecoli_pos_lo, int ecoli_pos_hi):
	"""
	Similar to gather_reads_BowTie but in preparation for subsampling,
	simply goes through the file, marking down all file locations of eligible reads

	Returns -- a list of (file locations, quality read len) of reads that can be used for subsampling
	"""
	cdef int i, len_read, phred, cur_tell
	eligible = []
	f = open(filename)
	while True:
		cur_tell = f.tell()
		line = f.readline()
		if not line: break
		name, strand, ref_seq_id, offset, read, quality, junk = line.strip().split('\t', 6)
		len_read = len(read)
		if refmap.gap_map[ref_seq_id][offset] < ecoli_pos_lo or\
				refmap.gap_map[ref_seq_id][offset] > ecoli_pos_hi:
					continue
		for i in range(len_read+1):
			if i == len_read:
				break
			phred = ord(quality[i]) - BOWTIE_PHRED_OFFSET
			assert 0 <= phred <= 40
			if phred < phred_cutoff:
				break
		if i >= min_length:
			eligible.append((cur_tell,i))
	f.close()
	return eligible

def subsample_reads_BowTie(object filename, object refmap, np.ndarray[np.int64_t, ndim=2] seqvec, object eligible, int size):
	"""
	Given eligible which is a list of (file pos, read len) from subsample_reads_BowTie_prepare
	Randomly sample <size> reads and store the nt counts in seqvec (which should be cleared before called!)
	"""
	import random
	cdef int N, i, j, pos, ind
	N = len(eligible)
	f = open(filename)
	for i in random.sample(xrange(N), min(N,size)):
		f.seek(eligible[i][0])
		name, strand, ref_seq_id, offset, read, quality, junk = f.readline().strip().split('\t', 6)
		offset = int(offset)
		for j in range(eligible[i][1]):
			pos = refmap.gap_map[ref_seq_id][offset+j]
			nt = read[j]
			if nt not in NT_MAPPING:
				nt = 'N'
			ind = NT_MAPPING[nt]
			seqvec[ind, pos] += 1
	f.close()

def gather_reads_inhouse(object filename, object refmap, np.ndarray[np.int64_t, ndim=2] seqvec, int phred_cutoff, int min_length, int max_degen, int ecoli_pos_lo, int ecoli_pos_hi):
	"""
	Currently UNUSED and UNTESTED
	max_degen -- max number of mismatches of the read to the refseq
	"""
	cdef int i, j, pos, ind, len_read, mismatch
	cdef int used=0, discard=0
	with open(filename) as f:
		aligned = cPickle.load(f)[1] # pickle is (unaligned, aligned), just look at aligned
	for m in aligned.itervalues(): # m is a BowTieMatch
		len_read = len(m.read)
		offset = int(m.offset)
		if refmap.gap_map[m.ref_seq_id][offset] < ecoli_pos_lo or\
				refmap.gap_map[m.ref_seq_id][offset] > ecoli_pos_hi:
#			raw_input("discarding read becuz mapped pos " + str(refmap.gap_map[ref_seq_id][offset]))
			discard += 1
			continue
		for i in range(len_read+1):
			if i == len_read:
				break
			if ord(m.quality[i]) - BOWTIE_PHRED_OFFSET < phred_cutoff:
				break
		if i >= min_length:
			j = 0
			mismatch = 0
			to_use = True
			while j < i:
				nt_a = m.read[j]
				nt_b = refmap.fasta[m.ref_seq_id].seq[m.offset+j]
				if (nt_a in ('T','U') and nt_b not in ('T','U')) or \
				   (nt_b in ('T','U') and nt_a not in ('T','U')) or \
				   (nt_a not in ('T','U') and nt_b not in ('T','U') and nt_a!=nt_b):
					mismatch += 1
				if mismatch > max_degen:
#					raw_input("not using becuz mismatch maxed:\n{0}\n{1}".format(\
#							refmap.fasta[m.ref_seq_id].seq[m.offset:(m.offset+j)],\
#							m.read[:j]))
					discard += 1
					to_use = False
					break
				j += 1
			if not to_use:
				continue
			used += 1
#			print("using the first {0} bases:\n(ref){1}\n(cur){2}".format(i,\
#					refmap.fasta[m.ref_seq_id].seq[m.offset:(m.offset+i)], m.read[:i]))
			for j in range(i):
				pos = refmap.gap_map[m.ref_seq_id][m.offset+j]
				nt = m.read[j]
				if nt not in NT_MAPPING:
					nt = 'N'
				ind = NT_MAPPING[nt]
				seqvec[ind, pos] += 1
		else:
#			print("discarding read {0} becuz useful read len {1}".format(m.read, i))
			discard += 1
	return used, discard
	
