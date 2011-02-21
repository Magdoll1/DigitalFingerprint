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
#SAMPLES = ['1418-1', '1418-4', '1419-1', '1419-4', '1420-1', '1420-4']

if os.popen("hostname").read().startswith('rhino'):
	print >> sys.stderr, "switching silo due to rhino...."
	os.environ['silo'] = "/shared/silo_researcher/Lampe_J/Gut_Bugs/"

from Solexa_settings import L2
valid_DI_pos = filter(lambda i: L2[i]%1==0 and 358<=L2[i]<=514, xrange(520))

filename = os.environ['silo'] + "/FH_Meredith/data/16S_090630/{0}/alignedSILVA100justHumanCrapV3curated1crap{0}.trimmedB2M30.fastq.out"
inhouse = os.environ['silo'] + "/FH_Meredith/data/16S_090630/{0}/unalignedSILVA100*inhouse_aligned.pickle*"
refmap = Read.RefMap(os.environ['PCODE'] + '/Silva/SILVA100_justHumanCrap.1crap_masked.V3region.fna.gap_map.bz2', 520, os.environ['PCODE'] + '/Silva/SILVA100_justHumanCrap.1crap_masked.V3region_ungapped.fna')
phred_cutoff = 10
min_length = 30
max_degen = 10
ecoli_lo = L2.index(340)
ecoli_hi = L2.index(484)

def main(file_iter, output_df_filename, log_f):
	log_f.write("phred cutoff:{0}\n".format(phred_cutoff))
	log_f.write("min length:{0}\n".format(min_length))
	log_f.write("max degen (if used):{0}\n".format(max_degen))
	log_f.write("use ecoli range {0}-{1}\n".format(ecoli_lo, ecoli_hi))
	f = open(output_df_filename, 'w')
	dfwriter = DF.DFWriter(f)
	for sample,file in file_iter:
		print >> sys.stderr, "processing {0}.........".format(sample)
		seqvec = np.zeros((5,520), dtype=np.int)
		# --------------- for in-house aligned ------------- #
		for file in glob.iglob(inhouse.format(sample)):
			print >> sys.stderr, "unzipping {0}.......".format(file)
#			if file.endswith('.bz2'):
#				os.system("bunzip2 " + file)
#				file = file[:-4]
			used, discarded = hello.gather_reads_inhouse(file, refmap, seqvec, phred_cutoff, min_length, max_degen, ecoli_lo, ecoli_hi)
			log_f.write("FILE:{0} USED:{1} DISCARDED:{2}\n".format(file, used, discarded))
#			os.system("bzip2 " + file)
		#  ---------------- for BowTie-aligned ---------------#
#		used, discarded = hello.gather_reads_BowTie(file, refmap, seqvec, phred_cutoff, min_length, ecoli_lo, ecoli_hi)
#		print >> sys.stderr, "used:", used, "discarded:", discarded
#		log_f.write("FILE:{0} USED:{1} DISCARDED:{2}\n".format(file, used, discarded))
		df = Read.ReadDF(sample, refmap)
		df.len = 520
		df.assign_vec(seqvec)
		dfwriter.write(df)
		runner = DiversityIndexRunner()
		di=runner.run(df, method='Simpson', threshold=0, vec_pre_normalized=False, ignoreN=True)[valid_DI_pos]
		print("{0},{1}".format(sample,",".join(map(str,di))))
	f.close()

def random_seed_test(d='/mnt/silo/silo_researcher/Lampe_J/Gut_Bugs/FH_Meredith/data/16S_090630/1411-1',\
		pattern='1411-1.trimmedB2M30.bowtie_random_seed_*.aligned_out',\
		output_df_prefix='1411-1.trimmedB2M30.bowtie_random_seed'):
	"""
	For all files under <d> with pattern *trimmedB2M30.bowtie_random_seed_*.aligned_out
	run it through main() to get the DF and DI
	"""
	file_iter = [(os.path.basename(f),f) for f in glob.iglob(d + '/' + pattern)]
	output_df_filename = output_df_prefix + ".phred{0}min{1}.DF".format(phred_cutoff, min_length)
	main(file_iter, os.path.join(d,output_df_filename))

if __name__ == "__main__":
	file_iter = [(s,filename.format(s)) for s in SAMPLES] # for 14-individual
	f = open('test.out_bowtie_phred10min30ecoli340to484.log', 'a+')
	main(file_iter, 'test.out_bowtie_phred10min30ecoli340to484.DF',f)
#	random_seed_test()
#	df_dict = {}
#	for df in DF.DFReader(open('test.DF')):
#		df_dict[df.name] = df
#	for df in DF.DFReader(open('test.out_bowtie_phred10min30.DF')):
#		df_dict[df.name].assign_vec(df_dict[df.name].vec+df.vec)
#	f = open('test.out_bowtie_inhouse_phred10min30maxdegen5_v2.DF','w')
#	w = DF.DFWriter(f)
#	w.writes(df_dict.itervalues())
#	f.close()
	
