import os,re,sys
import functools
from collections import defaultdict
from cPickle import *
from Bio import SeqIO

"""
After running RDP classifier online, the detail result can be downloaded and
will have the format:

Rfull_i21c1;  ; Root; 100%; Bacteria; 100%; Bacteroidetes; 100%; Bacteroidetes; 98%; Bacteroidales; 98%; Prevotellaceae; 97%; Prevotella; 65%

Here are scripts that will parse these details.

The 10 levels should be:
	0: Root
	1: Kingdom
	2: Phylum
	3: Class
	4: Sub-Class (-idae)
	5: Order (-ales)
	6: Sub-Order (-ineae)
	7: Family (-aceae)
	8: Genus (-monas, -ium, -ella, -bacter, -bacterium, -bacillus, -coccus, -spira, -vibrio)
	9: Species

"""
assembled_id_rex = re.compile('R\S+_i\d+c(\d+)l\d+')

RDP_LEVEL = ['domain', 'phylum', 'class', 'subclass', 'order', 'suborder', 'family', 'subfamily', 'genus', 'species']

GENUS_ENDINGS = ['monas', 'ium', 'ella', 'bacter', 'bacterium', 'bacillus', 'coccus', 'spira', 'vibrio', 'phagus']
GENUS_ENDINGS += ['Xylophilus', 'Blautia']

KNOWN_CLASS = ['Bacilli', 'Clostridia', 'Erysipelotrichi', 'Mollicutes', 'Verrucomicrobiae', 'Deinococci']
SAME_PHYLUM_CLASS_NAME = [\
			'Actinobacteria',\
			'Fusobacteria',\
			'Aquificae',\
			'Deferribacteres',\
			'Spirochaetes',\
			'Thermodesulfobacteria',\
			'Gemmatimonadetes',\
			'Chlamydiae',\
			'Chloroflexi',\
			'Thermotogae',\
			'Chrysiogenetes',\
			'Acidobacteria',\
			'Spirochaetes',\
			]

extra_rex = re.compile('(\S+) (\d+)')
nasty_rex = re.compile('([a-zA-Z]+(_\d+)*)')

def is_family(tax):
	return tax.endswith('aceae') or tax.startswith('Incertae')

def get_smart_tax_level(tax, last_lv):
	if last_lv < 2: return last_lv + 1
	elif last_lv == 2: # last was Phylum
		if tax.endswith('ia') or tax in KNOWN_CLASS or tax in SAME_PHYLUM_CLASS_NAME: return 3 # Class
		elif tax.endswith('idae'): return 4 # Sub-Class
		elif tax.endswith('ales'): return 5 # Order
		elif tax.endswith('ineae'): return 6 # Sub-Order
		elif is_family(tax): return 7 # Family
		elif any( [tax.endswith(x) for x in GENUS_ENDINGS] ): return 8 # Genus
		else: print >> sys.stderr, "ODD: should be Class but tax is {0}".format( tax )
	elif last_lv == 3: # last was Class
		if tax.endswith('idae'): return 4 # Sub-Class
		elif tax.endswith('ales'): return 5 # Order
		elif tax.endswith('ineae'): return 6 # Sub-Order
		elif is_family(tax): return 7 # Family
		elif any( [tax.endswith(x) for x in GENUS_ENDINGS] ): return 8 # Genus
		else:
			print >> sys.stderr, "ODD: should be Order but tax is {0}".format( tax )
	elif last_lv == 4: # last was Sub-Class
		if tax.endswith('ales'): return 5
		elif tax.endswith('ineae'): return 6 # Sub-Order
		elif is_family(tax): return 7 # Family
		elif any( [tax.endswith(x) for x in GENUS_ENDINGS] ): return 8 # Genus
		else:
			print >> sys.stderr, "ODD: should be Order but tax is {0}".format( tax )
	elif last_lv == 5: # last was Order
		if tax.endswith('ineae'): return 6 # Sub-Order
		elif is_family(tax): return 7 # Family
		elif tax.endswith('oideae'): return 8 # Sub-Family, pretend Genus?!
		elif any( [tax.endswith(x) for x in GENUS_ENDINGS] ): return 8 # Genus
		else:
			print >> sys.stderr, "ODD: should be Family but tax is {0}".format( tax )
	elif last_lv == 6: # last was Sub-Order
		if is_family(tax): return 7 # Family
		elif tax.endswith('oideae'): return 8 # Sub-Family, pretend Genus!?
		elif any( [tax.endswith(x) for x in GENUS_ENDINGS] ): return 8 # Genus
		else:
			print >> sys.stderr, "ODD: should be Family but tax is {0}".format( tax )
	elif last_lv == 7: # last was Family
		if tax.endswith('oideae'):
			print >> sys.stderr, "ODD: should be Genus but tax is {0}".format( tax )
		elif any( [tax.endswith(x) for x in GENUS_ENDINGS] ): return 8 # Genus
		elif tax.find(' sp.') > 0:
			print >> sys.stderr, "ODD: should be Genus but tax is {0}".format( tax )
		else: return 8
	elif last_lv == 8: # last was Genus
		if tax.find(' sp.') > 0:
			print >> sys.stderr, "ODD: should be Genus but tax is {0}".format( tax )
		else: return 9
	
	return last_lv + 1

def split_slv_tax(tax, splitter):
	"""
	Given an slv-style taxonomy (could be separated by either ;, / or _), return 
	"""
	parsed = []
	for t in tax.split(splitter):
		if t.count('_') == 0: # no chance of those nasty tax ex: Clostridia_1
			parsed.append( t )
		else:
			# for any tax like Clostridia_1, simply remove the _1....
			p = [r.group(1) for r in re.finditer(nasty_rex, t)]
			if len(p) == 1: # was not nasty
				i = t.rfind('_')
				if i >= 0 and str.isdigit(t[i+1:]):
					parsed.append( t[:i] )
				else:
					parsed.append( t )
			else:
				for x in p:
					i = x.rfind('_')
					if i >= 0 and str.isdigit(x[i+1:]): 
						parsed.append( x[:i] )
					else:
						parsed.append( x )
	return parsed

def parse_slv_tax(tax, splitter=';'):
	"""
	Given an slv-style taxonomy, ex: 
	    Bacteria/Actinobacteria/Acidimicrobidae_Acidimicrobiales_Acidimicrobineae/marine group

	Return a parsed hierarchy dict where key is tax level (0-9), value is the name of that tax
	   ex: {0: 'Root', 1: 'Bacteria', 2: 'Actinobacteria', 3: 'Actinobacteria', 
	        4: 'Acidimicrobidae', .....}
	"""
	result = {0: 'Root'}
	last_lv = 0

	# remove duplicates
	xx = split_slv_tax(tax, splitter )
	if len(xx) != len(set(xx)):
		i = 0
		while i < len(xx):
			if xx[i] in xx[:i]: xx.pop(i)
			else: i += 1
#		print >> sys.stderr, "DUP REMOVED, now", xx
	for t in xx:
		last_lv = get_smart_tax_level(t, last_lv)
		result[ last_lv ] = t

	return result

def iterate_tax_assignment(t):
	"""
	Given <t> which is in format:
		Root;1.;Bacteria;1.;Firmicutes;0.98.....

	iteratively return (level, assignment, confidence)
	"""
	tt      = t.split(';')
	last_lv = -1

	for i in xrange(0, len(tt)-1, +2):
		name    = tt[i].replace('\"', '') # strip the quotes
		conf    = float(tt[i+1])
		last_lv = get_smart_tax_level(name, last_lv)
		yield (last_lv, name, conf)

def parse_RDPdetail(detail_filename, abundance_lambda, min_confidence, get_detail=False):
	"""
	<detail_filename> should be the RDP classifier output in the form:
		
		EF698023\t\tRoot\tno rank\t1.0\tBacteria\tdomain\t1.0\t"Firmicutes"\tphylum\t1.0\t"Clostridia"...

	abundance_lambda -- given the seq ID (ex: Rfull_i21c1) returns the abundance of the seq
	min_confidence   -- minimum confidence bootstrap threshold, in floating point (ex: 0.8)

	returns a list T where T[level] is taxonomy --> count, level = 0,1,....9
	"""
	if get_detail:
		detail_info = {} # will be seq ID --> tax_level --> assignment
	else:
		tally_at_level = [ defaultdict(lambda: 0) for x in xrange(10)] # RDP can have up to 10 levels

	with open(detail_filename) as f:
		for line in f:
			line = line.strip().split('\t')
			id = line[0]
			if get_detail:
				detail_info[id] = {}

			#for lv,name,conf in iterate_tax_assignment(str(r.seq)):
			for i in xrange( (len(line)-5)/3 ):
				name, rank, conf = line[5+3*i], line[6+3*i], float(line[7+3*i])
				name = name.replace('"', '')
				conf = float(conf)
				try:
					lv   = RDP_LEVEL.index( rank )
				except ValueError:
					raise Exception, "unfound rank {0}!!!".format( rank )
		 		if conf >= min_confidence:
					if get_detail:
						detail_info[ id ][ lv ] = name
					elif lv < 10:
						tally_at_level[ lv ][ name ] += abundance_lambda(id)


	for i in xrange(10): tally_at_level[i] = dict(tally_at_level[i])

	if get_detail:
		return detail_info
	else:
		return tally_at_level

def print_tally_at_level(t, f, N):
	kk = t.keys()
	kk.sort()
	for k in kk:
		f.write("{0}\t{1}\n".format(k, t[k]))
	f.write("UNCLASSIFIED\t{0}\n".format( N - sum(t.itervalues()) ))

def print_tally(T, f):
	# ignores level 0 which is just Root
	for lv in xrange(1, 10):
		if len(T[lv]) > 0:
			f.write("---level {0}---\n".format(RDP_LEVEL[lv]))
			print_tally_at_level(T[lv], f, N=sum(T[0].itervalues()))
			f.write("\n")

def create_DBfile_for_RDPclassifier(nds_filename, output_filename):
	"""
	Given an input file that has the format per line:
		id \t tax_slv
	
	Write out to <output_filename> the tax file to be used for RDP classifier training
	with format:
		taxid*taxon name*parent taxid*depth*rank"
	"""
	tax_seen = set()
	tax_parsed = {}
	collect = defaultdict(lambda: set()) # depth --> set of seen taxons at that depth
	tax_to_id = {} # taxon name --> taxid

	with open(nds_filename) as handle:
		for line in handle:
			id, tax = line.strip().split('\t')
			if tax in tax_seen: continue
			tax_seen.add( tax )
			p = parse_slv_tax( tax )
			kk = p.keys()
			for k in kk:
				if k not in REAL_DEPTH: del p[k]
			for depth,name in p.iteritems(): collect[depth].add( name )
			tax_parsed[tax] = p

	id = 0
	for depth in xrange(10): # RDP or SLV should have up to 10 depth/levels
		for name in collect[depth]:
			tax_to_id[ name ] = id
			id += 1
	
	return tax_to_id, tax_parsed

def write_DBfile_for_RDPclassifier(output_filename, tax_parsed):
	import itertools
	# make tax_to_id
	pool_at_depth = defaultdict(lambda: set())
	for v in tax_parsed.itervalues():
		for depth, name in v.iteritems():
			pool_at_depth[depth].add( name )
	id = 0
	tax_to_id = defaultdict(lambda: {})
	for depth,names in pool_at_depth.iteritems():
		for name in names:
			tax_to_id[depth][name] = id
			id += 1

	name_written = defaultdict(lambda: set ())
	with open(output_filename, 'w') as f:
		for p in tax_parsed.itervalues():
			depths = p.keys()
			depths.sort(reverse=True)
			for i in xrange(len(depths)-1):
				depth  = depths[i]
				parent = tax_to_id[ depths[i+1] ][ p[depths[i+1]] ]
				name   = p[depth]
				if name in name_written[depth]: continue
				name_written[depth].add( name )
				
				f.write("{taxid}*{name}*{parent}*{depth}*{rank}\n".format(\
						taxid  = tax_to_id[ depth ][ name ],\
						name   = name,\
						parent = parent,\
						depth  = REAL_DEPTH[depth],\
						rank   = REAL_RANK[depth]))

def fixing_some_bad_tax_parsed(tax_parsed):
	"""
	To get RDP Classifier trainer to be happy with our tax_parsed,
	we have to make sure that every entry has all depths in REAL_DEPTH, 
	so here are some fixes:

	(1) if phylum is Actinobacteria and class is missing, make it Actinobacteria
	(2) append 'Unclassified' to tail depths
	"""
	SAME_PHYLUM_CLASS_NAME = [\
			'Actinobacteria',\
			'Fusobacteria',\
			'Aquificae',\
			'Deferribacteres',\
			'Spirochaetes',\
			'Thermodesulfobacteria',\
			'Gemmatimonadetes',\
			'Chlamydiae',\
			'Chloroflexi',\
			'Thermotogae',\
			'Chrysiogenetes',\
			'Acidobacteria',\
			'Spirochaetes',\
			]

	for v in tax_parsed.itervalues():
		if 8 in v and (v[8]=='Incertae Sedis' or v[8].startswith('uncultured') or v[8].startswith('Unclassified') or v[8].startswith('unlcassified') or v[8].startswith('environ')) :
			v[8] = 'Unclassified_' + v[7]
			print >> sys.stderr, "fixed to", v
			
		if set(v.keys()) == set(REAL_RANK_NEEDED): continue
		
		if 3 not in v and 2 in v:
			if v[2] in SAME_PHYLUM_CLASS_NAME:
				v[3] = v[2]

		if 5 in v and 3 not in v and v[5] in ('Thermales', 'Deinococcales'):
			v[3] = 'Deinococci'

		if 2 in v and v[2].startswith('Candidate division'):
			e = v[2][len('Candidate division'):].strip()
			v[3] = 'Unclassified_Class_' + e
			v[5] = 'Unclassified_Order_' + e
			v[7] = 'Unclassified_Family_' + e
			v[8] = 'Unclassified_Genus_' + e

		if 8 not in v and len(v.keys())==5:
			v[8]='Unclassified_Genus_' + v[7]

		if 3 in v and v[3].startswith('Subsection') and 8 in v:
			v[5]='Unclassified_Order_'+v[8]
			v[7]='Unclassified_Family_'+v[8]

		if 7 not in v and 8 not in v and len(v)==4:
			v[7] = 'Unclassified_Family_' + v[5]
			v[8] = 'Unclassified_Genus' + v[5]

		if 8 in v and 7 not in v and len(v)==5:
			v[7] = 'Unclassified_Family_' + v[5]

		print >> sys.stderr, "fixed to", v


def create_FASTAfile_for_RDPclassifier(nds_filename, fasta_filename, tax_map, output_filename):
	"""
	Given the NDS file (id --> tax), the fasta file, tax_map (raw tax --> parsed depths),
	write to <output_filename> the format for RDP Classifier training
	"""
	id_to_raw_tax = {}
	with open(nds_filename) as f:
		for line in f:
			id,tax = line.strip().split('\t')
			id_to_raw_tax[id] = tax
	
	f = open(output_filename, 'w')
	with open(fasta_filename) as handle:
		for r in SeqIO.parse(handle, 'fasta'):
			try:
				p = tax_map[ id_to_raw_tax[r.id] ]
			except:
				print >> sys.stderr, "ignoring for", r.id
				continue
			tax_str = ""
			for depth in xrange(10):
				if depth in p: tax_str += p[depth] + ';'
			f.write(">{0} {1}\n{2}\n".format(r.id, tax_str, r.seq))
	f.close()

def parse_assembled_id(id):
	m = assembled_id_rex.match( id )
	return int(m.group(1))

if __name__ == "__main__":
	file1 = 'data/RDPclassifier/20090908_Fei4414aligned_justV3.fna_download.txt'
	file2 = 'data/RDPclassifier/BowTie_processed.SILVA100justHumanCrapV3curated20090908_V3frags51_for_Richard.REASSEMBLED.fna_download.txt'

	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option("-f", "--file", help="fasta filename")
	parser.add_option("-w", "--weightedID", help="used weighted seq ID, default False", \
			action="store_true", dest="IDweighted", default=False)
	parser.add_option("-c", "--confidence", help="confidence level cutoff, default 0.10", \
			type=float, default=0.10)

	options, args = parser.parse_args()

	lamb = lambda x: int(x[x.find('c', x.find('_'))+1:])

	if options.IDweighted:
		T = parse_RDPdetail(options.file, parse_assembled_id, options.confidence)
	else:
		print >> sys.stderr, "ID NOT WEIGHTED"
		T = parse_RDPdetail(options.file, lambda x: 1, options.confidence)

	
	with open(options.file + '_parsed', 'w') as f:
		print_tally(T, f)
	
	with open(options.file + '_parsed.pickle', 'w') as f:
		dump(T, f)
