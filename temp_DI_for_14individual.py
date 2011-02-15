import os,re,sys,glob
import numpy as np
import hello
import Read
import DF
from DiversityIndex import DiversityIndexRunner

SAMPLES = ['1411-1', '1411-4', '1412-1', '1412-4', '1413-1', '1413-4', \
        '1414-1', '4414-1', '1414-4', '4414-4', '1415-1', '1415-4', \
        '1416-1', '1416-4', '1417-1', '1417-4', '1418-1', '1418-4', \
        '1419-1', '1419-4', '1420-1', '1420-4']

from Solexa_settings import L2
valid_DI_pos = filter(lambda i: L2[i]%1==0 and 358<=L2[i]<=514, xrange(520))

filename = "/mnt/silo/silo_researcher/Lampe_J/Gut_Bugs/FH_Meredith/data/16S_090630/{0}/alignedSILVA100justHumanCrapV3curated1crap{0}.trimmedB2M30.fastq.out"
inhouse = "/mnt/silo/silo_researcher/Lampe_J/Gut_Bugs/FH_Meredith/data/16S_090630/{0}/unalignedSILVA100*inhouse_aligned.pickle*"
refmap = Read.RefMap(os.environ['PCODE'] + '/Silva/SILVA100_justHumanCrap.1crap_masked.V3region.fna.gap_map.bz2', 520)
phred_cutoff = 20
min_length = 30

f = open('test.out_bowtie_phred20min30.DF', 'w')
dfwriter = DF.DFWriter(f)
for sample in SAMPLES:
	print >> sys.stderr, "processing {0}.........".format(sample)
	seqvec = np.zeros((5,520), dtype=np.int)
#	for file in glob.iglob(inhouse.format(sample)):
#		print >> sys.stderr, "unzipping {0}.......".format(file)
##		if file.endswith('.bz2'):
#			os.system("bunzip2 " + file)
#			file = file[:-4]
#		hello.gather_reads_inhouse(file, refmap, seqvec, phred_cutoff, min_length)
#		os.system("bzip2 " + file)
	g = hello.gather_reads_BowTie(filename.format(sample), refmap, seqvec, phred_cutoff, min_length)
	df = Read.ReadDF(sample,refmap)
	df.len = 520
	df.assign_vec(seqvec)
	dfwriter.write(df)
	f.flush()
	runner = DiversityIndexRunner()
	di=runner.run(df, method='Simpson', threshold=10)[valid_DI_pos]
	print("{0},{1}".format(sample,",".join(map(str,di))))
f.close()
