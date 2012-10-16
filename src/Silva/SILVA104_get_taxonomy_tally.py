import os, sys, csv
import parse_taxRDP
from collections import defaultdict

names = [line.strip().split()[0] for line in open('SILVA104.fece_augmented.fasta.gap_map')]

level_count = defaultdict(lambda: defaultdict(lambda: 0)) # level --> tax --> count
for r in csv.DictReader(open('SILVA104.nds'), delimiter='\t'):
	if r['name'] in names:
		for lv, tax in parse_taxRDP.parse_slv_tax(r['tax_slv']).iteritems():
			level_count[lv][tax] += 1

