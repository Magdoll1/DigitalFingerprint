import sys
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from utils import versatile_open
open = versatile_open.versatile_open

class FastaReader:
	"""
	This is meant to substitute for the Bio.SeqIO.to_dict method since some fasta files
	are too big to fit entirely to memory. The only requirement is that every id line
	begins with the symbol >. It is ok for the sequences to stretch multiple lines.
	The sequences, when read, are returned as Bio.SeqRecord objects.

	Example:
		r = FastaReader('output/test.fna')
		r['6C_49273_NC_008578/2259031-2259297'] ==> this shows the SeqRecord
	"""
	def __init__(self, fasta_filename, rename_dups=False):
		self.f = open(fasta_filename)
		self.d = {} # seq id --> file location
		self.dup_counter = 0 # we would NOT being using this unless rename_dups is True
		
		while 1:
			line = self.f.readline()
			if len(line) == 0: break
			if line.startswith('>'):
				id = line.strip()[1:]
				if id in self.d:
					if rename_dups:
						# we handle dups by RENAMING them
						# this is ok ONLY if we never have external references to the IDs!!
						self.dup_counter += 1
						new_id = id + '.' + str(self.dup_counter)
						print >> sys.stderr, "duplicate id {0} found. Renaming to {1}!!".format(id, new_id)
						id = new_id
					else:
						# we don't handle dups and this is NOT ok
						raise KeyError, "duplicate id {0}!!".format(id)

				self.d[id] = self.f.tell()

	def __iter__(self):
		for id in self.d:
			yield id

	def __getitem__(self, k):
		if k not in self.d:
			raise KeyError, "key {0} not in dictionary!".format(k)
		self.f.seek(self.d[k])
		content = ''
		for line in self.f:
			if line.startswith('>'):
				break
			content += line.strip().replace(' ','') # remove all blanks in case there are some
		return SeqRecord(Seq(content), id=k)

	def iterkeys(self):
		return self.__iter__()

	def keys(self):
		return self.d.keys()

if __name__ == "__main__":
	f = FastaReader(sys.argv[1])
