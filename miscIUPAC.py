import os,re,sys
from sets import ImmutableSet
import itertools

IUPAC_RNA_MAP = \
		{'N':ImmutableSet(['N']),\
		'A':ImmutableSet(['A']),\
		'U':ImmutableSet(['U']),\
		'G':ImmutableSet(['G']),\
		'C':ImmutableSet(['C']),\
		'R':ImmutableSet(['A','G']),\
		'Y':ImmutableSet(['C','U']),\
		'S':ImmutableSet(['G','C']),\
		'W':ImmutableSet(['A','U']),\
		'K':ImmutableSet(['G','U']),\
		'M':ImmutableSet(['A','C']),\
		'B':ImmutableSet(['C','G','U']),\
		'D':ImmutableSet(['A','G','U']),\
		'H':ImmutableSet(['A','C','U']),\
		'V':ImmutableSet(['A','C','G'])}

IUPAC_DNA_MAP = \
		{'N':ImmutableSet(['N']),\
		'A':ImmutableSet(['A']),\
		'T':ImmutableSet(['T']),\
		'G':ImmutableSet(['G']),\
		'C':ImmutableSet(['C']),\
		'R':ImmutableSet(['A','G']),\
		'Y':ImmutableSet(['C','T']),\
		'S':ImmutableSet(['G','C']),\
		'W':ImmutableSet(['A','T']),\
		'K':ImmutableSet(['G','T']),\
		'M':ImmutableSet(['A','C']),\
		'B':ImmutableSet(['C','G','T']),\
		'D':ImmutableSet(['A','G','T']),\
		'H':ImmutableSet(['A','C','T']),\
		'V':ImmutableSet(['A','C','G'])}

IUPAC_RNA_MAP_rev = dict( [(y,x) for (x,y) in IUPAC_RNA_MAP.items()] )
IUPAC_RNA_MAP_rev[ ImmutableSet(['A','C','U','G']) ] = 'N'
IUPAC_DNA_MAP_rev = dict( [(y,x) for (x,y) in IUPAC_DNA_MAP.items()] )
IUPAC_DNA_MAP_rev[ ImmutableSet(['A','C','T','G']) ] = 'N'

def get_IUPAC_RNA_code(nts):
	if 'N' in nts:
		nts.remove('N')
	if len(nts) == 0:
		return 'N'
	joined = reduce(lambda acc,nt: IUPAC_RNA_MAP[nt].union(acc), nts, ImmutableSet())
	return IUPAC_RNA_MAP_rev[ joined ]

def get_IUPAC_DNA_code(nts):
	if 'N' in nts:
		nts.remove('N')
	if len(nts) == 0:
		return 'N'
	joined = reduce(lambda acc,nt: IUPAC_DNA_MAP[nt].union(acc), nts, ImmutableSet())
	return IUPAC_DNA_MAP_rev[ joined ]

def can_pair(x, y):
	if x == 'T': x = 'U'
	if y == 'T': y = 'U'
	set_x = IUPAC_RNA_MAP[x]
	set_y = IUPAC_RNA_MAP[y]
	for a,b in itertools.product(set_x, set_y):
		if (a == 'A' and b == 'U') or \
		   (a == 'U' and b == 'A') or \
		   (a == 'C' and b == 'G') or \
		   (a == 'G' and b == 'C') or \
		   (a == 'G' and b == 'U') or \
		   (a == 'U' and b == 'G'):
			   return True
	return False

if __name__ == "__main__":
	print can_pair(sys.argv[1], sys.argv[2])
