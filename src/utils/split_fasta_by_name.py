import os
import sys
from Bio import SeqIO

def split_fasta_by_name(input_fasta, split_func, suffix='.fna', extra_criteria=lambda x: True):
	"""
	Given <input_fasta>, split it into different fasta files
	that are identified by <split_func> which is a mapping func of seq id --> name 
	
	ex: all TW9_xxxx seqs go into TW9.fna while all TW1_xxxx seqs go to TW1.fna
	"""
	handles = {}
	for r in SeqIO.parse(open(input_fasta), 'fasta'):
		name = split_func(r.id) 
		out = name + suffix
		if out not in handles:
			if os.path.exists(out):
				raise Exception, "file {0} already exists! Not OK! Abort!".format(out)
			else:
				f = open(out, 'w')
				handles[out] = f
		else:
			f = handles[out]
		if extra_criteria is None or extra_criteria(r):
			f.write(">{id}\n{seq}\n".format(id=r.id, seq=r.seq))
	
	for f in handles.itervalues():
		print >> sys.stderr, "output written to", f.name
		f.close()

if __name__ == "__main__":
	first_underline = lambda x: x[:x.find('_')]

	def split_func(x):
		i = x.rfind('.')
		if i > 0: return x[:i]
		else: return x[:x.rfind('_')]

	def extra_criteria(seq_record):
		return len(seq_record.seq) >= 200

	split_fasta_by_name(sys.argv[1], \
			split_func=first_underline, \
			suffix='.fna', \
			extra_criteria=None)
		
