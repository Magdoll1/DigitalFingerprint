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
inhouse_eligible = os.environ['silo'] + "/FH_Meredith/data/16S_090630/{0}/inhouse_eligible_phred{1}min{2}maxdegen{3}ecoli{4}to{5}.pickle"
refmap = Read.RefMap(os.environ['PCODE'] + '/Silva/SILVA100_justHumanCrap.1crap_masked.V3region.fna.gap_map.bz2', 520, os.environ['PCODE'] + '/Silva/SILVA100_justHumanCrap.1crap_masked.V3region_ungapped.fna')
phred_cutoff = 10
min_length = 30
max_degen = 3
ecoli_lo = L2.index(340) # this is used for removing reads that map to pos ...338, 339 of V3
ecoli_hi = L2.index(484) # this is used for removing reads that map to pos 484, 485, ...of V3

SUBSAMPLE_SIZES = [10, 20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480, 40960, 81920, 163840]

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
#			if file.endswith('.bz2'):
#				os.system("bunzip2 " + file)
#				file = file[:-4]
			used, discarded = hello.gather_reads_inhouse(file, refmap, seqvec, phred_cutoff, min_length, max_degen, ecoli_lo, ecoli_hi)
			print >> sys.stderr, file, used, discarded
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

def subsampling_prepare_inhouse(file_iter):
	from cPickle import *
	for sample,file in file_iter:
		eligible = []
		print >> sys.stderr, "processing {0}.........".format(sample)
		for file in glob.iglob(inhouse.format(sample)):
			eligible += hello.subsample_reads_inhouse_prepare(file, refmap, phred_cutoff, min_length, max_degen, ecoli_lo, ecoli_hi)
		print >> sys.stderr, "{0} total eligible in-house reads for sample {1}".format(len(eligible), sample)
		f = open(inhouse_eligible.format(sample, phred_cutoff, min_length, max_degen, L2[ecoli_lo], L2[ecoli_hi]), 'w')
		print >> sys.stderr, "picking reads to ", f.name
		dump(eligible, f)
		f.close()

def subsampling_BowTie_n_inhouse(file_iter, di_method):
	"""
	Run subsampling for BowTie+inhouse (inhouse must already have been pre-processed for eligibility)
	Simply prints out per line: <sample>,<size>,<comma-separated DI>
	"""
	from cPickle import *
	runner = DiversityIndexRunner()
	seqvec = np.zeros((5,520), dtype=np.int)
	for sample,file in file_iter:
		eligible_bowtie = hello.subsample_reads_BowTie_prepare(file, refmap, phred_cutoff, min_length, ecoli_lo, ecoli_hi)
		eligible_inhouse = load(open(inhouse_eligible.format(sample, phred_cutoff, min_length, max_degen, L2[ecoli_lo], L2[ecoli_hi])))
		p = len(eligible_bowtie)*1./(len(eligible_bowtie)+len(eligible_inhouse))
		for size in SUBSAMPLE_SIZES:
			print >> sys.stderr, sample, size
			seqvec[:] = 0
			hello.subsample_reads_BowTie(file, refmap, seqvec, eligible_bowtie, int(size*p))
			hello.subsample_reads_inhouse(refmap, seqvec, eligible_inhouse, phred_cutoff, min_length, size-int(size*p))
			df = Read.ReadDF(sample, refmap)
			df.len = 520
			df.assign_vec(seqvec)
			di=runner.run(df, method=di_method, threshold=0, vec_pre_normalized=False, ignoreN=True)[valid_DI_pos]
			print("{0},{1},{2}".format(sample,size,",".join(map(str,di))))

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

def jackknifing_tree(file_pattern, di_method):
	"""
	Given a pattern for the list of subsampled DI files
	(each file should have per-line format <sample>,<size>,<comma-separated DI>
	 and have one of the sizes be 'real')
	Run clustering and count the differences between subsampled and real trees
	Returns: tree (size-->sample-->list of trees), symmetric_difference (size-->list of diffs),
	         robinson_foulds_distance (size-->list of diffs)
	"""
	import numpy as np
	from clustering import Cluster
	import dendropy
	dTree = lambda x: dendropy.Tree.get_from_string(x, "newick")
	samples = None
	sizes = None
	trees = {}
	for file in glob.iglob(file_pattern):
		print >> sys.stderr, "reading subsamplined DI file {0}....".format(file)
		d = {}
		with open(file) as f:
			for line in f:
				sample, size, di = line.strip().split(',', 2)
				if size not in d:
					d[size] = {}
				d[size][sample] = np.array(map(float, di.split(',')))
		if sizes is None:
			sizes = d.keys()
			sizes.sort()
			samples = d[sizes[0]].keys()
			samples.sort()
		for size, di_dict in d.iteritems():
			c = Cluster(None)
			c.init_from_di_list(di_dict, method=di_method, threshold=0)
			c.run_till_end()
			try:
				trees[size].append(dTree(str(c.trees[0])))
			except KeyError:
				trees[size] = [dTree(str(c.trees[0]))]

	# tally (1) symmetric differences (edge weight ignored)
	# (2) robinson_foulds_distance (edge weight considered)
	# 'real' is the size that is the full pool that we compare all other trees to
	sym_diff = {}
	rob_diff = {}
	for size in sizes:
		if size == 'real': continue
		t_real = trees['real'][0]
		sym_diff[size] = [t_real.symmetric_difference(t) for t in trees[size]]
		rob_diff[size] = [t_real.robinson_foulds_distance(t) for t in trees[size]]

	return trees, sym_diff, rob_diff

if __name__ == "__main__":
	file_iter = [(s,filename.format(s)) for s in SAMPLES] # for 14-individual
	#subsampling_prepare_inhouse(file_iter)
	result = jackknifing_tree('useThisFlora/*v2_iter*', di_method='Entropy')
	sys.exit(-1)
	f = open('test.out_inhouse_phred10min30maxdegen5ecoli340to484.log', 'a+')
	main(file_iter, 'test.out_inhouse_phred10min30maxdegen5ecoli340to484.DF',f)
#	random_seed_test()
	df_dict = {}
	for df in DF.DFReader(open('test.out_inhouse_phred10min30maxdegen5ecoli340to484.DF')):
		df_dict[df.name] = df
	for df in DF.DFReader(open('test.out_bowtie_phred10min30ecoli340to484.DF')):
		df_dict[df.name].assign_vec(df_dict[df.name].vec+df.vec)
	f = open('test.out_bowtie_inhouse_phred10min30maxdegen5_v3.DF','w')
	w = DF.DFWriter(f)
	w.writes(df_dict.itervalues())
	f.close()
	
