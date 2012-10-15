import os, sys
from collections import defaultdict
from Bio import SeqIO
"""
A quick script to randomly sample at most N pyrosequences

I'm using this to generate a smaller test file
 from DAR_HCT_Bv3__ALL_DATASETS.Export.fa
 for demonstration

 the sample ids are like:
 >EZ0R7OU02IVY5I|DAR_HCT_Bv3|TD1|0.700%
 means it's from sample TD1
"""
SEQ_MAX = 200

def main(filename, output_prefix):
	count = defaultdict(lambda: 0)
	for r in SeqIO.parse(open(filename), 'fasta'):
		#id,ignore,sample,ignore = r.id.strip().split('|', 3)
		id, sample = r.id, r.id[2:6]
		if count[sample] < SEQ_MAX:
			os.system("echo \">{id}\" >> {out}_{sample}.fna".format(\
					id=id, out=output_prefix, sample=sample))
			os.system("echo \"{seq}\" >> {out}_{sample}.fna".format(\
					seq=r.seq, out=output_prefix, sample=sample))
			count[sample] += 1

if __name__ == "__main__":
	main('16Spyro_28F_all.fna', '16Spyro_28F_200seq/')

