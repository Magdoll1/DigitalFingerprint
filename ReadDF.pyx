import numpy as np
cimport numpy as np
from Solexa_settings import BOWTIE_PHRED_OFFSET, NT_MAPPING, BowTieMatch
import cPickle

def gather_reads_BowTie(object filename, object refmap, np.ndarray[np.int64_t, ndim=2] seqvec, int phred_cutoff, int min_length):
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

def gather_reads_inhouse(object filename, object refmap, np.ndarray[np.int64_t, ndim=2] seqvec, int phred_cutoff, int min_length):
	"""
	Currently UNUSED and UNTESTED
	"""
	cdef int i, j, pos, ind, len_read
	with open(filename) as f:
		aligned = cPickle.load(f)[1] # pickle is (unaligned, aligned), just look at aligned
	for m in aligned.itervalues(): # m is a BowTieMatch
		len_read = len(m.read)
		offset = int(m.offset)
		for i in range(len_read+1):
			if i == len_read:
				break
			if ord(m.quality[i]) - BOWTIE_PHRED_OFFSET < phred_cutoff:
				break
		if i >= min_length:
			#print("using the first {0} bases: {1}".format(i, read[:i]))
			for j in range(i):
				pos = refmap.gap_map[m.ref_seq_id][m.offset+j]
				nt = m.read[j]
				if nt not in NT_MAPPING:
					nt = 'N'
				ind = NT_MAPPING[nt]
				seqvec[ind, pos] += 1
		else:
			pass
	
