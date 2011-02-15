import os
import sys
from cPickle import *
import Read
from DF import DFWriter
from Solexa_settings import L2

if os.uname()[1].startswith('rhino'):
	print >> sys.stderr, "switching silo to /shared/silo_researcher/Lampe_J/Gut_Bugs/"
	os.environ['silo'] = '/shared/silo_researcher/Lampe_J/Gut_Bugs/'

"""
Currently stored in $silo/FH_Meredith/data/16S_090630 are folders 1411-1, 1411-4,....4414-4
 in each folder there should already be a inhouse_aligned.MERGED.pickle which is
 (aligned, read_count, phreds)

 where aligned is a list of (<read_id>, BowTieMatch object)
 and is the combo of all BowTie + in-house aligned stuff
"""
print >> sys.stderr, "reading refmap...."
refmap = Read.RefMap(os.environ['PCODE'] + '/Silva/SILVA100_justHumanCrap.1crap_masked.V3region.fna.gap_map.bz2', 520, os.environ['PCODE'] + '/Silva/SILVA100_justHumanCrap.1crap_masked.V3region_ungapped.fna')

SAMPLES = ['1411-1', '1411-4', '1412-1', '1412-4', '1413-1', '1413-4', \
		'1414-1', '4414-1', '1414-4', '4414-4', '1415-1', '1415-4', \
		'1416-1', '1416-4', '1417-1', '1417-4', '1418-1', '1418-4', \
		'1419-1', '1419-4', '1420-1', '1420-4']

DIR = os.path.join(os.environ['silo'], 'FH_Meredith/data/16S_090630')
PICKLE = 'inhouse_aligned.MERGED.pickle'
BOWTIE = "alignedSILVA100justHumanCrapV3curated1crap{sample}.trimmedB2M30.fastq.out"

SAMPLES = [sys.argv[1]]
output_df_filename = "14individual_Illumina_{0}.DF".format(sys.argv[1])

f = open(output_df_filename, 'w')
dfwriter = DFWriter(f)

for sample in SAMPLES:
	filename = os.path.join(DIR, sample, PICKLE)
	print >> sys.stderr, 'reading pickle', filename
	readdict = Read.ReadsDict(refmap)
	readdict.read_bowtie_match_pickle(filename)
	filename = os.path.join(DIR, sample, BOWTIE.format(sample=sample))
	print >> sys.stderr, 'reading bowtied', filename
	readdict.read_bowtie_output(filename)
	print >> sys.stderr, "remove reads mapped to pos 338, 339, >= 485"
	for i in xrange(4):
		del readdict.M[i]
	print >> sys.stderr, "total {0} reads".format(readdict.count)
	readdf = Read.ReadDF(name=sample, refmap=refmap)
	for read in readdict:
		readdf.add_read_to_vec_using_ref(read)#readdf.add_read_to_vec(read)
	dfwriter.write(readdf)
f.close()
