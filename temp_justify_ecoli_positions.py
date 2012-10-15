import os, re, sys
from Bio import SeqIO
from collections import defaultdict
from Solexa_settings import L2
from SILVA import Ecoli1542_SILVA100

seq_fasta_filename = '/shared/silo_researcher/Lampe_J/Gut_Bugs/FH_Meredith/data/DigFinger_wrapup_work/refDB/SILVA104.fece_augmented.fasta'
ref = SeqIO.read(open('Silva/ARB_ECOLI.fna'), 'fasta').seq.tostring().replace('U', 'T')

def constancy(ref, seq, i):
	"""
	simple 1-mer scenario: can we move ref[i] to i-1, ... or i+1...?
	without worsening the alignment quality?
	"""
	for j in xrange(i-1, -1, -1):
		if seq[j] != '-':
			break
		if ref[j] == seq[i]:
#			print >> sys.stderr, "possible altenative aln (ecoli {0}->{1}):".format(L2[i],L2[j])
#			print >> sys.stderr, ref[j:i+1]
#			print >> sys.stderr, seq[i]+seq[j+1:i]+'-'
			return True
	for j in xrange(i+1, len(ref)):
		if seq[j] != '-':
			break
		if ref[j] == seq[i]:
#			print >> sys.stderr, "possible alternative aln2(ecoli {0}->{1}):".format(L2[i], L2[j])
#			print >> sys.stderr, ref[i:j+1]
#			print >> sys.stderr, '-'+seq[i+1:j]+seq[i]
			return True
	return False

def main(input_filename, ref):
	"""
	For every ecoli and non-ecoli position, tally
	(1) number (%) of gaps
	(2) constancy of alignment
	"""
#	gap = defaultdict(lambda: 0)
#	cos = defaultdict(lambda: 0)
	nonecoli = defaultdict(lambda: 0) # r.id --> # of non-ecoli positions
	hasecoli = defaultdict(lambda: 0) # r.id --> # of ecoli positions
	for r in SeqIO.parse(open(input_filename), 'fasta'):
		seq = r.seq.tostring().replace('U', 'T').replace('.', '-')
		print >> sys.stderr, r.id
		for i, x in enumerate(seq):
			if x == '-': # is gap
				pass #gap[i] += 1
			else:
				if i in Ecoli1542_SILVA100:
					hasecoli[r.id] += 1
				else:
					nonecoli[r.id] += 1
				#cos[i] += constancy(ref, seq, i)
#		raw_input("press to continue")
	return hasecoli, nonecoli #return gap, cos

print >> sys.stderr, "convert all U to T!!!"
y, n = main(seq_fasta_filename, ref)
