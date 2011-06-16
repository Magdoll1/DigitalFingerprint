import os, sys, itertools

IUPAC_RNA_MAP = \
		{'N':frozenset(['N']),\
		'A':frozenset(['A']),\
		'U':frozenset(['U']),\
		'G':frozenset(['G']),\
		'C':frozenset(['C']),\
		'R':frozenset(['A','G']),\
		'Y':frozenset(['C','U']),\
		'S':frozenset(['G','C']),\
		'W':frozenset(['A','U']),\
		'K':frozenset(['G','U']),\
		'M':frozenset(['A','C']),\
		'B':frozenset(['C','G','U']),\
		'D':frozenset(['A','G','U']),\
		'H':frozenset(['A','C','U']),\
		'V':frozenset(['A','C','G'])}

IUPAC_DNA_MAP = \
		{'N':frozenset(['N']),\
		'A':frozenset(['A']),\
		'T':frozenset(['T']),\
		'G':frozenset(['G']),\
		'C':frozenset(['C']),\
		'R':frozenset(['A','G']),\
		'Y':frozenset(['C','T']),\
		'S':frozenset(['G','C']),\
		'W':frozenset(['A','T']),\
		'K':frozenset(['G','T']),\
		'M':frozenset(['A','C']),\
		'B':frozenset(['C','G','T']),\
		'D':frozenset(['A','G','T']),\
		'H':frozenset(['A','C','T']),\
		'V':frozenset(['A','C','G'])}

IUPAC_RNA_MAP_rev = dict( [(y,x) for (x,y) in IUPAC_RNA_MAP.items()] )
IUPAC_RNA_MAP_rev[ frozenset(['A','C','U','G']) ] = 'N'
IUPAC_DNA_MAP_rev = dict( [(y,x) for (x,y) in IUPAC_DNA_MAP.items()] )
IUPAC_DNA_MAP_rev[ frozenset(['A','C','T','G']) ] = 'N'

def get_IUPAC_RNA_code(nts):
	if 'N' in nts:
		nts.remove('N')
	if len(nts) == 0:
		return 'N'
	joined = reduce(lambda acc,nt: IUPAC_RNA_MAP[nt].union(acc), nts, frozenset())
	return IUPAC_RNA_MAP_rev[ joined ]

def get_IUPAC_DNA_code(nts):
	if 'N' in nts:
		nts.remove('N')
	if len(nts) == 0:
		return 'N'
	joined = reduce(lambda acc,nt: IUPAC_DNA_MAP[nt].union(acc), nts, frozenset())
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
