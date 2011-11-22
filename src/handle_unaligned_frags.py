from Solexa_settings import *
from difflib import SequenceMatcher as SM, Match
import miscParses
from hybrid_structure import HybridCons

refDBdir  = '/home/etseng/FH_Meredith/SolexaReads/refDBs/'
refDBname = 'SILVA100_justHumanCrap.1crap_masked.V3region_ungapped_nonredundant.Nmasked.fna'

d = {}
with open(os.path.join(refDBdir,refDBname)) as handle:
	for r in SeqIO.parse( handle, 'fasta' ):
		d[r.id] = SM(None, r.seq.data)

HC = HybridCons()

BestFragMatch = namedtuple('BestFragMatch', 'ref_seq_id i bpq match_count match_phred strand')
parseBLAST    = miscParses.parseNCBIBLASTline

IMMEDIATE_BEST_MATCH_COUNT = 47
MINIMUM_MATCH_COUNT        = 35 # if the best match is below this, consider it unaligned

def score_match_by_basequality(s1, s2, q):
	"""
	Given two matching sequences (which must be of same length)
	and q, a Solexa quality string (Phred = ord()-64), return
	a matching score weighted by the base qualities
	"""
	assert len(s1) == len(s2)
	score = 0
	for i in xrange(len(s1)):
		if s1[i]==s2[i]: score += ord(q[i]) - SOLEXA_PHRED_OFFSET
	return score

def blast_frags(frags_filename, seed_length, top_strand_only=False):
	"""
	Run BLAST on {frags_filename} given {seed_length}

	Returns query_name --> [ list of NCBI-blast matches ]

	SIDE EFFECT: deletes both {frags_filename} and the BLAST output
	"""
	cmd  = "blastall -p blastn -v 99999 -b 99999 -G -999999 -E -999999 -m 8 -a 4 -v 20 -b 20 "
	cmd += "-W {0} ".format( seed_length )
	cmd += "-d {0} ".format( refDBname )
	cmd += "-i {0} ".format( frags_filename )

	# make sure the file is not there first
	if os.path.exists( frags_filename + '.out' ):
		os.remove( frags_filename + '.out' )
	
	# we do it once on the + strand
	# and then once on the - strand
	os.system( "touch " + frags_filename )
	#print >> sys.stderr, "running blast {0} on + strand".format( frags_filename )
	os.system( cmd + " -S 1 >> {0}.out".format( frags_filename ) )
	if not top_strand_only:
		#print >> sys.stderr, "running blast {0} on - strand".format( frags_filename )
		os.system( cmd + " -S 2 >> {0}.out".format( frags_filename ) )

	# parse the BLAST output, remember the positions here are all 1-based
	result = defaultdict(lambda: [])
	for line in open( frags_filename + '.out' ):
		p = parseBLAST(line)
		result[p.id1].append( p )
	# delete the files as we don't need it anymore
	os.remove( frags_filename )
	os.remove( frags_filename + '.out' )

	return result

def find_match_for_frag(frag_id, plus_frag, frag_quality, blasted, aligned, unaligned):
	"""
	Given the {blasted} results for {frag}, go through each hit, extending it 
	appropriately to see what the match_count (total number of identical overlap) is,
	and find the best refseq match for {frag}.

	If the best refseq match count >= MINIMUM_MATCH_COUNT, then the match info
	is put into {aligned} as frag_id --> BowTieMatch

	Otherwise it is appended to {unaligned} as (frag_id, frag, frag_quality)

	As a side-effect, converts all qualities to BowTie-style
	"""
	minus_frag = Seq(plus_frag).reverse_complement().data

	best = BestFragMatch(None, None, bpq=0., match_count=0, match_phred=0, strand=None)

#	print >> sys.stderr, "for frag {0} there are {1} cands".format(plus_frag, len(blasted))

	for pm in blasted:
		if pm.gaps > 0:
			#print >> sys.stderr, "ignoring gaps for", pm
			continue # we don't handle gaps now
		assert pm.start1 < pm.end1
		size = pm.end1 - pm.start1 + 1
		match_count = size - int(pm.pcmismatch)

		if pm.start2 < pm.end2:
			frag        = plus_frag
			match       = Match( a=pm.start2-1, b=pm.start1-1, size=size )
			strand      = '+'
		else: # reverse-complement the read
			frag        = minus_frag
			frag_quality= frag_quality[::-1]
			match       = Match( a=pm.end2-1, b=len(frag)-pm.end1, size=size )
			strand      = '-'

		if match.a < match.b:
			#unaligned.append('BEFORE')
			#return 'BEFORE'
			continue # ABORT CASE 1: the read aligns BEFORE the beginning of ref seq

		# just faking a SequenceMatcher so we can re-use code below
		id = pm.id2 # ref_seq_id
		sm = d[id]
		sm.set_seq2(frag)
		match_phred = score_match_by_basequality(sm.a[match.a:match.a+size],\
				                                 frag[match.b:match.b+size],\
												 frag_quality[match.b:match.b+size])

		abort = False
		j = 1
		while match.b - j >= 0:
			if sm.a[ match.a-j ] == frag[ match.b-j ]:
				match_phred += ord(frag_quality[ match.b-j ]) - SOLEXA_PHRED_OFFSET
				match_count += 1
			j += 1
		j = 0
		while match.b + match.size + j < len(frag):
			if match.a + match.size + j >= len(sm.a):
				abort = True
				break # ABORT CASE 2: the read aligns BEYOND the end of the ref seq
			if sm.a[ match.a+match.size+j ] == frag[ match.b+match.size+j ]:
				match_phred += ord(frag_quality[ match.b-j ]) - SOLEXA_PHRED_OFFSET
				match_count += 1
			j += 1

		if abort:
			#unaligned.append('AFTER')
			#return 'AFTER'
			#print >> sys.stderr, "abort because aligns before"
			continue

		i = match.a - match.b
		q = HC.qualify_with_secondary_structure(frag, ref_seq_id=id, i=i, match_size=match.size)
#		else: (i < 0) # we've ABORTED this case
#			q = HC.qualify_with_secondary_structure(frag, ref_seq_id=id, i=0, match_size=match.size+i)

		if match_count >= IMMEDIATE_BEST_MATCH_COUNT and q == 1:
			best = BestFragMatch(id, i, q, match_count, match_phred, strand)
			break
		elif match_count >= MINIMUM_MATCH_COUNT and match_phred >= best.match_phred and q > best.bpq:
#			print >> sys.stderr, "{0} is a better than {1}".format(id, best)
			best = BestFragMatch(id, i, q, match_count, match_phred, strand)


	# as a side-effect, we're converting the Solexa-style quality string to BowTie-style!
	frag_quality = convert_phred_Solexa_to_BowTie(frag_quality)

	if best.ref_seq_id is None:
#		print >> sys.stderr, "{0} is still unaligned!".format(frag_id)
		unaligned.append( (frag_id,plus_frag,frag_quality) )
	else:
		if best.strand == '+':
			aligned[frag_id] = BowTieMatch(frag_id, best.strand, best.ref_seq_id, best.i, Seq(plus_frag), frag_quality, best.bpq)
		else:
			aligned[frag_id] = BowTieMatch(frag_id, best.strand, best.ref_seq_id, best.i, Seq(minus_frag), frag_quality, best.bpq)
		#print >> sys.stderr, "best for {0} is {1}".format(frag_id, aligned[frag_id])

def handle_unaligned(unaligned_filename, k, output_pickle=None):
	"""
	Given a file of unaligned reads (FASTQ format), try to align
	them by first BLAST-ing them (in bulk) against the refDB then use that as
	the candidate list for finding the best matching ref seq

	BLAST hits that indicate a match count of IMMEDIATE_BEST_MATCH_COUNT will
	 be immediately declared the best refseq match

	otherwise we go through the whole list of BLAST hits and output the best
	 one that has at least MINIMUM_MATCH_COUNT match count

	Returns
		unaligned --- a list of [(frag_name,frag,quality)]
		aligned   --- dict of frag_name --> BowTieMatch

	p.s. as a side-effect, all quality strings will be converted to BowTie-style
	"""
	BULK_SIZE = 100 # run the BLAST in bulks of <BULK_SIZE> reads

	unaligned, aligned = [], {}

	tmp_filename = os.tempnam()
	tmp_count = 0

	frag_quality = {}
	frag_seq     = {}

	with open(unaligned_filename) as handle:
		handle_EOF = os.stat( unaligned_filename ).st_size
		frag_f = open(tmp_filename, 'w')
		while handle.tell() != handle_EOF:
			id  = handle.readline().strip()[1:] # ex: @SOLEXA-1:....
			seq = handle.readline().strip()
			handle.readline() # same id again
			quality = handle.readline().strip()

			frag_f.write(">{0}\n{1}\n".format(id, seq))
			frag_quality[id] = quality
			frag_seq[id]     = seq
			tmp_count += 1

			if tmp_count >= BULK_SIZE:
				frag_f.close()

				result = blast_frags(frag_f.name, k)
				for frag_name,blasted in result.iteritems():
					find_match_for_frag(frag_name,\
							frag_seq[frag_name],\
							frag_quality[frag_name],\
							blasted, aligned, unaligned)

				tmp_count = 0
				frag_f = open(tmp_filename, 'w')
				frag_seq = {}
				frag_quality = {}
			
				# ignore, only pickle at the very end
				#if output_pickle is not None:
				#	with open(output_pickle, 'w') as pf: dump((unaligned, aligned), pf)

	# finish up the remainder
	if not frag_f.closed and frag_f.tell() > 0:
		frag_f.close()
		for frag_name,blasted in blast_frags(frag_f.name, k).iteritems():
				find_match_for_frag(frag_name,\
						frag_seq[frag_name],\
						frag_quality[frag_name],\
						blasted, aligned, unaligned)

	if output_pickle is not None:
		with open(output_pickle, 'w') as pf: dump((unaligned, aligned), pf)

	return unaligned, aligned

def compress_identical_reads(pickle_pattern):
	"""
	For all the pickles that conform to the "pickle pattern" ex: Solexa*inhouse_aligned.pickle
	read through each of them, collate the aligned reads, and when an identical read is encountered,
	increment the read_count

	Returns matches, read_count
	  where matches    -- dict of frag_name --> BowTieMatch
	        read_count -- read seq --> count
			phreds     -- dict of read seq --> np.array of accumulated phred scores
			unaligned  -- the sum of collection of unaligneds
	"""
	matches = {}
	read_count = {}
	phreds = {}
	unaligned = []

	dir,pattern = os.path.split( pickle_pattern )

	for file in fnmatch.filter( os.listdir(dir), pattern ):
		print >> sys.stderr, "reading pickle {0}....".format( file )
		with open( os.path.join(dir, file ) ) as handle:
			yet, aligned = load(handle)
			unaligned += yet
			for m in aligned.itervalues():
				q = [ ord(x)-BOWTIE_PHRED_OFFSET for x in m.quality ]
				if m.read.data in read_count:
					read_count[m.read.data] += 1
					phreds[m.read.data] += q
				else:
					matches[m.name] = m
					read_count[m.read.data] = 1
					phreds[m.read.data] = np.array( q ) 

	return matches, read_count, phreds, unaligned

if __name__ == "__main__":
	SEED_LENGTH = 4

	file = sys.argv[1]
	out  = file + '.inhouse_aligned.pickle'

	yet,done = handle_unaligned(file, SEED_LENGTH, out)
	with open(out, 'w') as f: dump((yet,done), f)

