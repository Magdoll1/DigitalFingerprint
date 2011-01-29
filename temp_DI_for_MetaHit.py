import os
import sys
import glob
import Read
from DF import DFWriter

from utils import versatile_open
open = versatile_open.versatile_open

"""
Harcoded temporary script file for collecting all the .bowtied files in Qin_MetaHit
(currently only the Spanish group: OC.XXX or V1.XXXx)
and output them into a DF file (designate output filename in command line)
"""

pattern = '/mnt/silo/silo_researcher/Lampe_J/Gut_Bugs/FH_Meredith/data/Qin_MetaHit/aligned/NotUsingForNow_MH/*.fq.bowtied'
refmap = Read.RefMap(os.environ['PCODE'] + '/Silva/SILVA104.fece_augmented.fasta.gap_map.bz2', aln_length=50000)

output_df_filename = "Qin_MetaHit.justDanish_dateseparate.DF"

f = open(output_df_filename, 'w')
dfwriter = DFWriter(f)

for file in glob.iglob(pattern):
	# file name ex: O2.UC-1_090112.fq.bowtied
	name = os.path.basename(file)[:-len('.fq.bowtied')]
	if len(os.popen("grep \"{0}\" {1}".format(name, output_df_filename)).read().strip()) > 0:
		print >> sys.stderr, "{0} already done, ignore!".format(name)
		continue
	print >> sys.stderr, "reading {0} for DF writing....".format(name)
	readdict = Read.ReadsDict(refmap)
	readdict.read_bowtie_output(file)
	readdf = Read.ReadDF(name, refmap)
	for read in readdict:
		readdf.add_read_to_vec(read)
	dfwriter.write(readdf)

f.close()
