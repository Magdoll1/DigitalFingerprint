from Solexa_settings import *
import miscIUPAC
import bisect
import handle_unaligned_frags as Shelp
import random

class Read:
	"""
	The unique identifier (id) will be the very first lane_tile_x_y given
	(which will happen when the object is first initialized)

	lane_tile_x_y_s --- stores a list of lane:tile:x:y, which are different
	                    reads that have the same sequence
						ex: ['1:1:35:12', '1:1:134:158', '1:4:24:100'....]

	copy --- the number of "available" copies. this is used during the
	         reassembly process so DO NOT USE THIS to get how many identical
			 sequence reads there were from the original export file!!!
			 instead, use len(lane_tile_x_y_s), since lane_tile_x_y_s
			 is NEVER altered after the initial export file set-up

	position --- the position this read mapped to on a particular reference
	             sequence. when the reference sequence changes, remember
				 to change this position accordingly
	
	end --- calculated as position + len(read) + 1, so it's +1 after actual end

	read --- a MutableSeq object, so if you want to compare two Read obj's sequences,
	         don't do r1.read==r2.read, use r1.read.data==r2.read.data, as
			 the former is comparing the two Seq objects, not the sequences 

	gapped_read --- a 0crap's gapped version of the read, by default is not
	         made (None), need to call make_gapped_read(ref_seq_id)

	gapped_end_1 --- same as above, default None
	ref_seq_id --- also default None, will be set when make_gapped_read is called
	"""
	def __init__(self, lane_tile_x_y, position, read, phreds):
		self.id = lane_tile_x_y
		self.lane_tile_x_y_s = [lane_tile_x_y]
		self.position = position
		self.end = position + len(read) # remember the end is +1 after the actual end
		self.read = read
		self.gapped_read = None
		self.gapped_end_1 = None
		self.ref_seq_id  = None
		self.copy = 1
		self.phreds = phreds # an np.array of quality score sums
	
		self.handle_weird_position()

	def __str__(self):
		return "[lane:tile:x:y]: {0}\n".format(self.lane_tile_x_y_s) + \
				"position: {0}\n".format(self.position) + \
				"copy: {0}\n".format(self.copy) + \
				"read: {0}\n".format(self.read)

	def __repr__(self):
		return "p{pos}c{copy}:{id}".format(pos=self.position,\
				copy=self.copy,id=self.id)

	def handle_weird_position(self):
		"""
		Currently, if self.position < 0,
		then trim it!! (TODO: is this the right thing todo?)
		"""
		if self.position < 0:
			delta = -self.position
			self.read = Seq(self.read.data[delta:])
			self.position = 0

	def make_gapped_read(self, ref_seq_id=None):
		"""
		Given the ref seq this read originally matched to,
		create 0crap's gapped version
		"""
		if ref_seq_id is None:
			if self.ref_seq_id is None:
				raise Exception, "to call make_gapped_read(), ref_seq_id must be given or already set!"
			else:
				ref_seq_id = self.ref_seq_id

		mapping = Shelp.HC.mapping[ref_seq_id]
		ungapped_i = mapping.index( self.position )
		ungapped_j = ungapped_i + len(self.read) - 1
		
		if ungapped_j < len(mapping):
			gapped_j = mapping[ ungapped_j ]
		else:
			print >> sys.stderr, "the read {0} goes beyond the refseq {1}, \
					trim short ({2})".format(ungapped_j, len(mapping), self.read.data)
			delta = ungapped_j - (len(mapping)-1)
			self.read = Seq(self.read.data[:-delta])
			assert mapping[ ungapped_i + len(self.read) - 1 ] == mapping[ -1 ]
			gapped_j = mapping[ -1 ]

		self.gapped_end_1 = gapped_j + 1
		self.gapped_read = Shelp.HC.refseqs[ ref_seq_id ].seq[ self.position : self.gapped_end_1 ]
		self.ref_seq_id = ref_seq_id

	def change_position(self, new_position):
		"""
		Important to change end if position is different!!
		"""
		self.position = new_position
		self.end = new_position + len(self.read)

	def add_copy(self, lane_tile_x_y):
		self.lane_tile_x_y_s.append( lane_tile_x_y )
		self.copy += 1

	def calc_identical_overlap(self, other):
		"""
		Return # of identical nucleotides at overlapping positions
		(note: positions where one or both have an 'N' is considered identical)

		If they don't overlap, will return 0
		"""
		score = 0
		for i in xrange(max(self.position, other.position),	min(self.end, other.end)):
			if self.read[i-self.position] == other.read[i-other.position] or \
				self.read[i-self.position] == 'N' or\
				other.read[i-other.position] == 'N':
				score += 1
		return score

class SeqVector:
	mapping = {'N':0, 'A':1, 'U':2, 'G':3, 'C':4}
	rev_mapping = ['N', 'A', 'U', 'G', 'C']
	gap_markers = ('.', '-')

	def __init__(self, start, end_1, position):
		"""
		SeqVector is used to gradually build a reassembly representation,
		once SeqVector is built, call consensus_seq() to get the reassembled sequence

		p.s. the sequence will be represented in RNA, not DNA

		start --- now defunct, just set to 0

		end_1 --- this is now the same as len

		position --- should be the position of the first read, as this is
		             used later to quickly identify potential mergings

		len --- the length of the reassembled sequence (NOT the number of reads!)

		vec --- a 5 x len numpy array where the entry (i, j) denotes the
		        frequency of nucleotide rev_mapping[i] at position j

		conss --- the consensus sequence, is not actually created until the
		          first time consensus_seq() is called

		conss_dirty_bit --- dirty bit flag for consensus_seq()

		read_count --- number of reads used to build this SeqVector, is automatically
		               incremented every time add_read or add_seq is called
		"""
		self.start = start
		self.end_1   = end_1  # +1 after the real end
		self.position= position
		self.len   = end_1 - start
		self.vec   = np.zeros((5, self.len), dtype=np.float)
		self.conss = ''
		self.conss_dirty_bit = True
		self.read_count = 0
		self.ids_at_position = []
	
	def add_read(self, read_obj, start_at_offset, copy=1):
		"""
		Add a Read obj at <offset> with abundance <copy>
		"""
		self.ids_at_position.append( (start_at_offset, read_obj.id) )
		self.read_count += copy
		rna_seq = read_obj.read.data.replace('T','U')
		for offset, x in enumerate( rna_seq ):
			if x in SeqVector.gap_markers:
				continue
			nt_s = miscIUPAC.IUPAC_RNA_MAP[ x ]
			w = copy*1. / len(nt_s)
			for nt in nt_s:
				try:
					self.vec[SeqVector.mapping[nt], start_at_offset+offset] += w
				except IndexError:
					print >> sys.stderr, "{0} with start_at_offset {1} offset {2} \
							exceeds length of {3}".format(read_obj.id, start_at_offset, offset, self.end_1)
					print >> sys.stderr, read_obj.gapped_read

		self.conss_dirty_bit = True # remember to set the dirty bit

	def add_seq(self, seq, start_at_offset, copy=1):
		"""
		Same effect as add_read, except that it's a seq
		"""
		self.read_count += copy
		for i,x in enumerate(seq):
			if x in SeqVector.gap_markers:
				continue
			nt_s = miscIUPAC.IUPAC_RNA_MAP[ x ]
			w = copy*1. / len(nt_s)
			for nt in nt_s:
				try:
					self.vec[SeqVector.mapping[nt], start_at_offset+i] += w
				except IndexError:
					pass
		self.conss_dirty_bit = True

	def remove_seq(self, seq, start_at_offset, copy=1):
		self.read_count -= copy
		for i,x in enumerate(seq):
			if x in SeqVector.gap_markers:
				continue
			nt_s = miscIUPAC.IUPAC_RNA_MAP[ x ]
			w = copy*1. / len(nt_s)
			for nt in nt_s:
				try:
					self.vec[SeqVector.mapping[nt], start_at_offset+i] -= w
				except IndexError:
					pass
		self.conss_dirty_bit = True

	def consensus_seq(self, gapped=False, percent_cutoff=.7):
		"""
		At each column, use code that represents >={percent_cutoff}% of the counts
		"""
		if not self.conss_dirty_bit:
			return self.conss

		seq = ''
		for col in xrange(self.len):
			total_count = self.vec[:, col].sum()
			nz = self.vec[:, col].nonzero()[0]
			if len(nz) == 0 and gapped:
				seq += '-'
			for size in xrange(1,len(nz)+1):
				done = False
				for st in itertools.combinations(nz, size):
					if sum( self.vec[x, col] for x in st ) >= total_count*percent_cutoff:
						seq += miscIUPAC.get_IUPAC_RNA_code( [SeqVector.rev_mapping[i] for i in st] )
						done = True
						break
				if done:
					break
		self.conss = seq
		self.conss_dirty_bit = False
		return seq

def remove_low_quality_for_matched(matches, read_count, phreds, min_phred_score, ditched_f):
	"""
	Remove any matches (and it's entries from {matches}, {read_count} and {phreds})
	that have 1 or more positions with a superPhred score < {min_phred_score}

	Returns 
	  count         -- total number of reads removed
	  count_unique  -- total number of unique reads removed
	"""
	count = count_unique = 0
	kk = matches.keys()
	for k in kk:
		m = matches[k]
		if any( x < min_phred_score for x in phreds[m.read.data] ):
			count += read_count[m.read.data]
			count_unique += 1
			ditched_f.write("@{id}\n{seq}\n+{id}\n{qual}\n".format( id=k, seq=m.read, \
					qual=m.quality ))
			del matches[k]
			del read_count[m.read.data]
			del phreds[m.read.data]
	return count, count_unique


def gather_reads_BowTie(export_filename, bested=None, read_count=None, phreds=None):
	"""
	Read <export_filename> which is the BowTie output file having fields:
	0. read name     1. strand     2. ref name
	3. 0-based offset (w.r.t + strand)     4.seq     5.quality string
	6. # of reserved (ignore in most cases)
	7. mismatch descriptors

	See http://bowtie-bio.sourceforge.net/manual.shtml#algn_out for output format.

	Returns a dict of read_name --> BowTieMatch tuple
	and a dict of read_seq --> count(copy)
	"""
	if bested is None: bested = {}
	if read_count is None: read_count = {}
	if phreds is None: phreds = {} # unique seq --> array of accumulated phred scores
	best_so_far = {}

	with open(export_filename) as handle:
		for line in handle:
			line = line.replace('\x00', '') # this seems to occassionally show up and screw things...?
			name,strand,ref_seq_id,offset,read,quality,junk = line.strip().split('\t',6)
			if name in bested:
				continue

			if read in read_count:
				read_count[read] += 1
				phreds[read] += [ ord(x)-BOWTIE_PHRED_OFFSET for x in quality ]
				continue
			else:
				read_count[read] = 1
				phreds[read] = np.zeros( len(read) )
				phreds[read] += [ ord(x)-BOWTIE_PHRED_OFFSET for x in quality ]

			read   = Seq(read)
			offset = int(offset)
			bpq = Shelp.HC.qualify_with_secondary_structure(read, ref_seq_id, offset, match_size=len(read))

			m = BowTieMatch(name, strand, ref_seq_id, offset, read, quality, bpq)

			# we don't need to worry about - strand. BowTie reverse-complements the reads for us!
			
			if bpq == 1.:
				bested[m.name] = m
				if m.name in best_so_far:
					del best_so_far[m.name]
			elif m.name in best_so_far:
				if bpq > best_so_far[m.name].bpq:
					best_so_far[m.name] = m
			else:
				best_so_far[m.name] = m

	bested.update( best_so_far )
	return bested, read_count, phreds

def makeMbyPos(matches, read_count, phreds):
	"""
	Returns MbyPos (or just abbreviated M) which is a dict of
	  (gapped) position --> list of Reads
	"""
	M = defaultdict(lambda: []) # (gapped) position --> list of Reads

	for m in matches.itervalues():
		gapped_i = Shelp.HC.mapping[ m.ref_seq_id ][ m.offset ]
		read_obj = Read(m.name, gapped_i, m.read, phreds[m.read.data])
		read_obj.ref_seq_id = m.ref_seq_id #read_obj.make_gapped_read( m.ref_seq_id )
		read_obj.copy = read_count[m.read.data]
			
		M[gapped_i].append( read_obj )
	return M

def makeUnalignedM(unaligned):
	"""
	Returns a mocked up version of MbyPos where
	 M[0] is (gapped) position --> list of Reads
	 the reads will have a position of inf since it's unknown
	"""
	read_count = {}
	phreds = {}
	matches = []
	M = {0: []}

	for id,seq,quality in unaligned:
		if seq in read_count:
			read_count[seq] += 1
			phreds[seq] += [ ord(x)-BOWTIE_PHRED_OFFSET for x in quality ]
		else:
			matches.append( (id,seq) )
			read_count[seq] = 1
			phreds[seq] = np.zeros( len(seq) )
			phreds[seq] += [ ord(x)-BOWTIE_PHRED_OFFSET for x in quality ]

	for id,seq in matches:
		read_obj = Read(id, float("inf"), Seq(seq), phreds[seq])
		read_obj.copy = read_count[seq]
		M[0].append( read_obj )

	del matches, read_count, phreds
	return M

def add_edge_thru_unaligned(M, unalignedM, G, hash_size, max_into, max_mismatch_func, tempHC=None):
	"""
	unalignedM should be a faked M where unalignedM[0] = list of Reads of unaligned
	(so position is inf)

	If exists at least one node u\in M s.t. u's seq is immediatel before it,
	make this unaligned a node v and add it in M (position is u's next pos)
	while removing it from unalignedM

	This WILL NOT WORK with trimmed reads of length < hash_size
	Also it only looks at ONE strand at a time, so remember to reverse-complement it
	 and run this again!
	"""
	def mix_in_place(s, t, l):
		if s.count('N') > 0 or t.count('N') > 0:
			new_s = list(s)
			new_t = list(t)
			for i in xrange(l):
				if s[i] == 'N' and t[i] != 'N': new_s[i] = t[i]
				elif t[i] == 'N' and s[i] != 'N': new_t[i] = s[i]
			return "".join(new_s), "".join(new_t)
		else: return s,t

	def count_mismatch(s, t, ps, pt, _from, _to_1):
		count = 0
		for i in xrange(_from, _to_1):
			if ps[i]>=BADPHRED_CUTOFF and\
			   pt[i]>=BADPHRED_CUTOFF and\
			   (s[i]!=t[i] and s[i]!='N' and t[i]!='N'): count += 1
		return count

	if tempHC is None: tempHC = Shelp.HC

	print("Remember you must run this twice, 2nd time with unalignedM contents reversed!")

	hash = defaultdict(lambda: []) # {hash_size}-key --> list of unaligned Reads having this prefix
	for x in unalignedM[0]:
		for p in permutate_N_in_seq( x.read.data[:hash_size], 0): hash[ p ].append( x )
	print >> sys.stderr, "we made a hash with {0} entries...".format( len(hash) )

	positions = M.keys()
	positions.sort()
	while len(positions) > 0:
		print >> sys.stderr, "doing position {0}".format(positions[0])
		V = M[ positions.pop(0) ]
		for v in V:
			for into in xrange(1, max_into+1):
				for p in permutate_N_in_seq( v.read.data[ into : hash_size+into ], 0):
					if p != v.read.data[ into : hash_size+into ]: continue
					# the hash we look at is hash[ p ]
					i = 0
					while i < len(hash[p]):
						x = hash[p][i]
						l = min(len(x.read.data), len(v.read.data))
						m = count_mismatch(v.read.data[into:l], x.read.data[:l-into], \
								v.phreds[into:l], x.phreds[into:l], hash_size, l-into)
						if max_mismatch_func( min(x.copy,v.copy), m ):
							# v --> x , so let's put x in M[next pos of v]
							#print >> sys.stderr, "found predecessor for unaligned {0} with into {1}".format(x.read.data, into)
							x_pos = v.position
							for junk in xrange(into):
								x_pos = tempHC.next_nongap_pos(v.ref_seq_id, x_pos)
							x.ref_seq_id = v.ref_seq_id
							x.position   = x_pos
							if M.has_key(x_pos):
								M[x_pos].append( x )
							else:
								M[x_pos] = [ x ]
								bisect.insort_left(positions, x_pos)
							# remove x from both unalignedM
							#unalignedM[0].remove( x ) # we'll remove it last
							# remove x (and all its variants if it contains 'N' from hash)
							for p2 in permutate_N_in_seq( x.read.data[:hash_size], 0):
								hash[ p2 ].remove( x )
							i -= 1 # rewind the index
							# in case of 'N'-containing, renew both seqs
							x.read.data, v.read.data = mix_in_place(x.read.data, v.read.data, l)
							# add the directed edge v --> x with weight being the matching length
							G.add_edge( v , x, (l-into-m, into) )
						i += 1

	# important unalignedM[0] removal here instead!
	# basically we just need to remove any entry with not-None-ref_seq_id
	i = 0
	while i < len(unalignedM[0]):
		if unalignedM[0][i].ref_seq_id is None: i += 1
		else: unalignedM[0].pop(i)

#	db = BowTie_writeout(unalignedM, xrange(0,1), 'head')
#	os.system("bowtie-build -q -r -o 1 {0} {0}_bowtie".format(db))
#	os.system("rm {0}_bowtie.rev.*".format(db)) # we don't want - hits
#	os.system("touch {0}_bowtie.rev.1.ebwt".format(db))
#	os.system("touch {0}_bowtie.rev.2.ebwt".format(db))
#
#	qr = BowTie_writeout(MbyPos, xrange(positions[0],positions[-1]+1), 'tail')
#
#	tmpout = os.tempnam()
#	os.system("bowtie -p 4 -n 0 -v 0 -a {0}_bowtie -f {1} {2}".format(db, qr, tmpout))
#	print >> sys.stderr, "db {0} query {1} tmpout {2}".format(db, qr, tmpout)
#
#	with open(tmpout) as f:
#		for line in f:
#			raw = line.strip().split()
#			u,strand,match_pos,v_name = name_to_read[raw[0]], raw[1], int(raw[3]), raw[2]
#			# UNFORTUNATELY BECAUSE WE TRIM READS THERE COULD BE BAD CASES HERE
#			# ex: ATT matches to TATTG, so match pos MUST BE 0 to indicate an
#			# immediate adjacent match!
#			# another possibility is that the unaligned matched to itself since we
#			# don't remove it from the BowTie db even though it's now a read
#			# another possibility is a reverse-strand match which is NOT okay here
#			# if we want to use unaligned's reverse we have to call this function again
#			if match_pos != 0 or strand!='+':
#				print >> sys.stderr, "bad case of matching {0} to {1}, ignore".format(u.id, v_name)
#				continue
#
#			v = name_to_read[ v_name ]
#			if v.ref_seq_id is None: # this means it's still an unaligned
#				# remove v from unalignedM
#				unalignedM[0].remove( v )
#				# we'll just pretend v is u's immediate next pos
#				v_pos = tempHC.next_nongap_pos(u.ref_seq_id, u.position)
#				v.position = v_pos
#				v.ref_seq_id = u.ref_seq_id #v.make_gapped_read( u.ref_seq_id )
#				if not MbyPos.has_key( v_pos ):
#					print >> sys.stderr, "inserting additional position {0}".format(v_pos)
#					MbyPos[ v_pos ] = []
#				MbyPos[ v_pos ].append( v )
#			else:
#				# It is entirely possible...that v was added prior to this position!
#				# which will make v.position > u.position but I guess it's okay....
#				print >> sys.stderr, "{0} has already been aligned, do nothing".format( v_name )
#
#	os.system("rm {0}*".format(qr)) 
#	os.system("rm {0}*".format(db)) 
#
#	print >> sys.stderr, "unalignedM is down to {0}".format( len(unalignedM[0]) )
#	print >> sys.stderr, "TIME ELAPSED: {0} sec".format(time.time()-start_t)

def allowed_mismatch_func(R, M, copy, mismatch):
	"""
	cutoffs should be a increasing-order list R denoting copy number
	 and a list M denoting allowed mismatches

	where the max # of mismatches allowed for a read with copy <= R[i] is M[i]
	"""
	return mismatch <= M[ bisect.bisect_left(R, copy) ]

def add_edge_thru_BowTieing(M, G, max_into, hash_size, delta, max_mismatch_func, tempHC=None):
	"""
	Uses BowTie to help with adding immediate edges u --> v 
	where u.read.data[i:] == v.read.data[:-i] where  1 <= i <= max_into

	"""
	def count_mismatch(s, t, ps, pt, _from, _to_1):
		#print >> sys.stderr, "comparing\n{0}\n{1}".format(s,t)
		count = 0
		for i in xrange(_from, _to_1):
			if ps[i]>=BADPHRED_CUTOFF and \
			   pt[i]>=BADPHRED_CUTOFF and \
			   (s[i]!=t[i] and s[i]!='N' and t[i]!='N'): count += 1
		#print >> sys.stderr, "mismatch: {0}".format(count)
		return count

	if tempHC is None: tempHC = Shelp.HC

	positions = M.keys()
	positions.sort()

	edge_count  = 0

	for pos in positions[:-1]:
		if pos == -1:
			max_next = positions[-1]
		else:
			max_next = pos + 1
			for ref_seq_id in tempHC.mapping:
				tmp = pos
				for junk in xrange(max_into):
					tmp = tempHC.next_nongap_pos(ref_seq_id, tmp)
				max_next = max( max_next, tmp )
		print >> sys.stderr, "max next of {0} is {1}".format(pos, max_next)

		hash = defaultdict(lambda: [])
		for x in M[pos]:
			for into in xrange(1, max_into+1):
				for p in permutate_N_in_seq( x.read.data[into:into+hash_size], 0 ):
					hash[ p ].append( (x, into) )
		print >> sys.stderr, "we made a hash of size {0}".format( len(hash) )

		tmp = [-1] + range( max(0, pos-delta), min(positions[-1]+1, max_next+delta+1) )
		for next in tmp:
			if next not in M: continue # sometimes M is no longer a defaultdict here...
			for v in M[next]:
				l_v = len( v.read.data )
				for p in permutate_N_in_seq( v.read.data[ : hash_size ], 0 ):
					for u, pos in hash[ p ]:
						#if G.has_edge(u, v): continue # already there, ignore
						l_u = len( u.read.data )
						# let's see if u --> v
						y = u.copy if u.copy < v.copy else v.copy
						if l_u < l_v+pos:
							m = count_mismatch(u.read.data[pos:],\
											   v.read.data[:l_u-pos],\
											   u.phreds[pos:],\
											   v.phreds[:l_u-pos],\
											   hash_size, l_u-pos)
							if max_mismatch_func( y, m ):
								G.add_edge( u, v, (l_u-pos-m, pos) )
#								raw_input("1:is a match!\n{0}\n{1} weight:{2}".format(u.read.data, ' '*(pos)+v.read.data,l_u-pos-m))
						elif l_u >= l_v+pos:
							m = count_mismatch(u.read.data[pos:pos+l_v],\
											   v.read.data,\
											   u.phreds[pos:pos+l_v],\
											   v.phreds,\
											   hash_size, l_v)
							if max_mismatch_func( y, m ):
								G.add_edge( u, v, (l_v-m, pos) )

		# do it especially again for all M[-1] reads that are not yet in G
		# (i.e. doesn't have an edge yet), but reverse-complemented
		print >> sys.stderr, "repeat for reverse-complement of unused M[-1]s...."
		for v in M[-1]:
			if v in G: continue # already in G, skip it
			v.read = v.read.reverse_complement()
			l_v = len( v.read.data )
			for p in permutate_N_in_seq( v.read.data[ : hash_size ], 0 ):
				for u, pos in hash[ p ]:
					#if G.has_edge(u, v): continue # already there, ignore
					l_u = len( u.read.data )
					# let's see if u --> v
					y = u.copy if u.copy < v.copy else v.copy
					if l_u < l_v+pos:
						m = count_mismatch(u.read.data[pos:],\
										   v.read.data[:l_u-pos],\
										   u.phreds[pos:],\
										   v.phreds[:l_u-pos],\
										   hash_size, l_u-pos)
						if max_mismatch_func( y, m ):
							G.add_edge( u, v, (l_u-pos-m, pos) )
							edge_count += 1
					elif l_u >= l_v+pos:
						m = count_mismatch(u.read.data[pos:pos+l_v],\
										   v.read.data,\
										   u.phreds[pos:pos+l_v],\
										   v.phreds,\
										   hash_size, l_v)
						if max_mismatch_func( y, m ):
							G.add_edge( u, v, (l_v-m, pos) )
							edge_count += 1
		print >> sys.stderr, "added {0} edges for reverse-complement".format(edge_count)					
		# flip it back
		for v in M[-1]:
			if v in G: continue # already in G, skip it
			v.read = v.read.reverse_complement()

	return edge_count

#		db = BowTie_writeout(M, xrange(pos, pos+1), 'tail', None)
#		os.system("bowtie-build -q -r -o 1 {0} {0}_bowtie".format(db))
#		print >> sys.stderr, "db for pos {0} is {1}".format(pos, db)
#
#		for next in xrange(min_next, max_next+delta+1):
#			print >> sys.stderr, "looking at the head of pos {0}".format(next)
#			for v in M[next]:
#				for line in os.popen("bowtie --quiet --norc -p 4 -n 0 -v 0 -a {0}_bowtie -c {1}".format(db, v.read.data[:hash_size-max_into+1])):
#					raw = line.strip().split()
#					strand, u, pos = raw[1], name_to_read[raw[2]], int(raw[3])
#					if strand != '+' or G.has_edge(u,v): continue
#					if pos < max_into:



#		db = BowTie_writeout(M, xrange(pos, pos+1), 'tail', None)
#		os.system("bowtie-build -q -r -o 1 {0} {0}_bowtie".format(db))
#		print >> sys.stderr, "db for pos {0} is {1}".format(pos, db)
#
#		for next in xrange(min_next, max_next+delta+1):
#			print >> sys.stderr, "looking at the head of pos {0}".format(next)
#			for v in M[next]:
#				for line in os.popen("bowtie --quiet --norc -p 4 -n 0 -v 0 -a {0}_bowtie -c {1}".format(db, v.read.data[:hash_size-max_into+1])):
#					raw = line.strip().split()
#					strand, u, pos = raw[1], name_to_read[raw[2]], int(raw[3])
#					if strand != '+' or G.has_edge(u,v): continue
#					if pos < max_into:
#						# we only hashed the first {hash_size} nucleotides so need to see if really match
#						l_u = len(u.read.data)
#						l_v = len(v.read.data)
#						l = min(l_u, l_v)
#						if l_u <= l_v and u.read.data[1+pos:l-pos] == v.read.data[:l-pos-1]:
#							print("is a match!\n{0}\n{1}".format(u.read.data, ' '*(1+pos)+str(v.read.data)))
#						elif l_u > l_v and u.read.data[1+pos:l+1-pos] == v.read.data:
#							print("is a match!\n{0}\n{1}".format(u.read.data, ' '*(1+pos)+str(v.read.data)))

def permutate_N_in_seq(seq, after):
	Ncount = seq.count('N', after)
	if Ncount == 0: return [seq]
	else:
		acc = []
		for nt in ('A','T','C','G'):
			new_seq = seq.replace('N', nt, 1)
			acc += permutate_N_in_seq(new_seq, seq.find('N', after))
		return acc

def BowTie_writeout(MbyPos, pos_range, head_or_tail, max_size):
	tmp_filename = os.tempnam()
	with open(tmp_filename, 'w') as f:
		for i in pos_range:
			if i not in MbyPos: continue
			for match in MbyPos[i]:
				if head_or_tail == 'head':
					if max_size is None: seq = match.read.data[:-1]
					else:                seq = match.read.data[:max_size]
					for w in permutate_N_in_seq(seq, 0):
						f.write(">{0}\n{1}\n".format(match.id, w))
				elif head_or_tail == 'tail':
					if max_size is None: seq = match.read.data[1:]
					else:                seq = match.read.data[ len(match.read.data)-max_size: ]
					for w in permutate_N_in_seq(seq, 0):
						f.write(">{0}\n{1}\n".format(match.id, w))
				elif head_or_tail == 'neither':
					for w in permutate_N_in_seq( str(match.read), 0):
						f.write(">{0}\n{1}\n".format(match.id, w))
				else:
					raise Exception, "{0} is invalid option for head_or_tail".format(head_or_tail)
	return tmp_filename

def make_result_from_out_filename(out_filename, topOrder):
	"""
	Given the output file produced by ShedSkin-ed version of iterativeMaximumPath
	which for each path found outputs:
	  L
	  offs
	  copy
	  path
	  //
	return the list of [(SeqVector, copy)....]
	"""
	result = []
	eof = os.stat(out_filename).st_size
	with open(out_filename) as f:
		while f.tell() != eof:
			L    = eval( f.readline() )
			offs = eval( f.readline() )
			copy = eval( f.readline() )
			path = eval( f.readline() )
			assert f.readline().strip() == '//'

			sv = SeqVector(0, L, topOrder[path[0]].position)
			for i,o in itertools.izip(path, offs):
				sv.add_read(topOrder[i], o, 1)

			result.append( (sv, copy) )
	return result

def entropy_index(result):
	"""
	Similar to diversity_index except that diversity_index uses Simpson Index
	and this uses the entropy information which is 
		sum_{s\in A,T,C,G} P(s)*log( P(s)/Q(s) )
	where Q is the background nucleotide distribution and
	      P is the current column's nucleotide distribution
	Right now we use log 2.
	"""
	log2 = lambda x: math.log(x, 2)
	# first we calculate the background distribution
	end = max(sv.end_1 for sv,copy in result)
	V = np.zeros((5, end))
	for sv, copy in result:
		V[:,sv.start:sv.end_1] += sv.vec * copy

	totalN = sum(sum(V[:,i] for i in xrange(end))) * 1.
	Q = [ sum(V[j,:]) / totalN for j in xrange(5) ]
	
	# now calculate per-column entropy
	D = {}
	for i in xrange(end):
		N = sum(V[:,i]) * 1.
		D[i] = 0
		for j in xrange(5):
			if V[j,i] != 0: 
				D[i] += (V[j,i] / N) * log2( (V[j,i] / N) / Q[j] )

	return D

def diversity_index(result, N_cutoff):
	"""
	<result> is a list of [(SeqVector, copy)....]
	returns D --- which is a dict of position --> Simpson Index
	  and   Ns --- a dict of position --> number of nucleotides used in calculation
	"""
	end = max( [sv.end_1 for sv,copy in result] )
	V = np.zeros((5,end))
	for sv,copy in result:
		V[:, sv.start:sv.end_1] += sv.vec * copy

	D = {}
	Ns = {}
	for i in xrange(end):
		N = sum(V[:, i])
		if N < N_cutoff:
			print >> sys.stderr, "skipping position {0} in DI-calculation becuz count {1} too low".format(i, N)
			continue
		Ns[i] = N
		if N == 0:
			D[i] = 0
		elif N == 1:
			D[i] = 0
		else:
			dd = 1-sum([s*(s-1) for s in V[:,i]])*1./(N*(N-1))
			D[i] = dd
	
	return D, Ns

def ARBaligned_into_SeqVector_list(fasta_filename, abundance_lambda):
	"""
	Read an ARB-aligned fasta file (where IDs are ex: R1420-1_AY974857i1c1s9e124)
	and make it into a list of [(SeqVector, copy)....]
	so we can use it later to calculate diversity index
	"""
	result = []

	with open(fasta_filename) as handle:
		for r in SeqIO.parse(handle, 'fasta'):
			copy = abundance_lambda(r.id)
			sv = SeqVector(0, len(r.seq), 0)
			ss = r.seq.tostring().replace('T','U')
			sv.add_seq( ss, 0 )
			result.append( (sv, copy) )

	return result

def ARBaligned_into_depth_of_coverage(abundance_lambda=lambda x:1, fasta_pattern='R*-*.done.ARB_aligned.fna',continuous=True):
	"""
	For each file that fits the fasta_pattern
	(likely R*-*.done.ARB_aligned.fna)
	get coverage, which is done by this:
	(1) for each seq, find it's first none-gap position (s)
	    and the last none-gap position (e)
	(2) s--e is covered once by this seq
	"""
	counts  = None
	for file in fnmatch.filter(os.listdir('.'), fasta_pattern):
		print >> sys.stderr, "processing {0}....".format(file)
		with open(file) as handle:
			for r in SeqIO.parse(handle, 'fasta'):
				if counts is None:
					# initialize counts
					counts = np.zeros( len(r.seq), dtype=np.int)
				sseq = r.seq.tostring()
				for i in xrange(len(sseq)):
					if sseq[i] not in ('-','.'):
						start = i
						break
				for i in xrange(len(sseq)-1, -1, -1):
					if sseq[i] not in ('-','.'):
						end = i
						break
				if continuous:
					counts[start:(end+1)] += abundance_lambda(r.id)
				else:
					counts[start] += abundance_lambda(r.id)
	return counts			

def depth_of_coverage(m):
	positions = m.keys()

	start = min(positions)
	end   = max(positions) + 51 # to be safe

	DoC = np.zeros(end)

	for p in positions:
		for read_obj in m[p]:
			DoC[read_obj.position] += read_obj.copy

	return DoC

def DoC_for_BowTieMatches(matches, read_count, L, continuous=True):
	"""
	Given <matches> which is a dict of read_name --> BowTieMatches
	Use the proper 0crap's ungapped-to-gapped mapping to get
	depth of coverage.

	NOTE: assumes we're doing this for V3 region
	"""
	DoC = np.zeros( L )

	for m in matches.itervalues():
		if m.offset < 0:
			i = Shelp.HC.mapping[m.ref_seq_id][0]
			print >> sys.stderr, "trimming head to {0}".format(i)
		else:
			i = Shelp.HC.mapping[m.ref_seq_id][m.offset]
		# be careful: mapping[offset+read_length-1]+1 not necessary eq mapping[offset+read_length]
		u_j_1 = m.offset + len(m.read)-1
		if u_j_1 >= len(Shelp.HC.mapping[m.ref_seq_id]):
			j_1 = Shelp.HC.mapping[m.ref_seq_id][-1] + 1
			print >> sys.stderr, "trimming tail to {0}".format(j_1)
		else:
			j_1 = Shelp.HC.mapping[m.ref_seq_id][m.offset + len(m.read)-1] + 1

		try:
			abundance = read_count[m.read.data]
		except:
			abundance = read_count[m.read]
		if continuous:
			DoC[ i : j_1 ] += abundance
		else:
			DoC[ i ] += abundance

	return DoC

def combine_dup_seq_vector(result):
	"""
	Given a list of [(SeqVector, copy)] see if exists two SeqVectors
	that are the same (i.e. same consensus_seq)

	This is possible only after we've combined incompletes....
	"""
	d = {} # consensus_seq --> (SeqVec, copy)
	for sv, copy in result:
		x = sv.consensus_seq()
		if x in d:
			d[x] = (d[x][0], d[x][1] + copy)
		else:
			d[x] = (sv, copy)
	return d.values()

def output_result(result, name, output_f, gapped=False):
	"""
	Given a list of [(SeqVector, copy)], output in Fasta format

	seq ID will be R{name}_i{index}c{copy}l{len}
	if copy number is fraction, take the nearest integer (but round to 1 if 0)
	"""
	for index, (sv,copy) in enumerate(result):
		seq = sv.consensus_seq(gapped)
		id = "R{name}_i{index}c{copy}l{len}".format(name=name, copy=max(1,int(round(copy))), index=index, len=len(seq))

		output_f.write(">{0}\n{1}\n".format( id, seq ))

def printDoC(DoC, L2, filename):
	"""
	Print two versions:
	[ NONZERO   ] print all entries that have a non-zero value
    [ SPACED    ] using L2, sum-up .5 positions, and print values for 338, 338.5, 339...etc
	"""
	X = DoC.nonzero()[0].tolist()
	f = open(filename, 'w')
	f.write( "NONZERO\n" )
	f.write( ",".join( [ str(L2[x]) for x in X ] ) + '\n' )
	f.write( ",".join( [ str(DoC[x]) for x in X ] ) + '\n' )

	sL2_str = ''
	sDoC    = []
	for i,x in enumerate(X):
		if i == 0:
			sDoC.append( DoC[x] )
			sL2_str += str(L2[x])
		elif L2[x] == L2[X[i-1]]: # same as last non-zero position, which means it's a .5
			sDoC[-1] += DoC[x]
		else:
			sL2_str += ',' + str(L2[x])
			sDoC.append( DoC[x] )

	f.write( "\nSPACED\n" )
	f.write( sL2_str + '\n' )
	f.write( ",".join( map(str, sDoC) ) + '\n' )
	f.close()

if __name__ == "__main__":
	assembled_id_rex = re.compile('R\S+_i\d+c(\d+)l\d+')
	def parse_assembled_id(id):
		m = assembled_id_rex.match( id )
		return int(m.group(1))

	for file in os.listdir('.'):
		if file.endswith('.fasta'):
			print >> sys.stderr, "calculating DI for {0}....".format(file)
			result = ARBaligned_into_SeqVector_list(file, lambda x:1)
			DI, Ns = diversity_index(result, 0)

			k = DI.keys()
			i = 0
			while i < len(k):
				if L2[k[i]]%1. != 0: k.pop(i)
				else: i += 1
			k.sort()
	
			print file
			print( ",".join( [ str(L2[i]) for i in k ] )  )
			print( ",".join( [ str(DI[i]) for i in k ] )  )
			print( ",".join( [ str(Ns[i]) for i in k ] )  )

