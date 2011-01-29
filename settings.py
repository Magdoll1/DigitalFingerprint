import os
import sys
import glob
import time
import datetime
from optparse import OptionParser
from Bio import SeqIO

def read_config(filename='config'):
	global SETTINGS
	SETTINGS = {}
	with open(filename) as f:
		for line in f:
			line = line.strip()
			if line.startswith('#') or len(line) == 0:
				continue
			key, val = line.strip().split()
			if key in SETTINGS:
				SETTINGS[key].append(val)
			else:
				SETTINGS[key] = [val]

	if 'REF' not in SETTINGS:
		SETTINGS['REF'] = os.path.join(os.getcwd(), 'refDB')
	else:
		SETTINGS['REF'] = SETTINGS['REF'][0]
	print >> sys.stderr, "reference DB is set to", SETTINGS['REF']
	if not os.path.exists(SETTINGS['REF']):
		print >> sys.stderr, "creating refDB directory since it doesn't exist"
		os.mkdir(SETTINGS['REF'])

	if 'GAP' not in SETTINGS:
		print >> sys.stderr, "no gap symbols!"
	else:
		print >> sys.stderr, "gap symbols are", SETTINGS['GAP']
	if 'LOG' not in SETTINGS:
		SETTINGS['LOG'] = 'log'
	else:
		SETTINGS['LOG'] = SETTINGS['LOG'][0]
	print >> sys.stderr, "log filename is set to", SETTINGS['LOG']
	if 'BOWTIE' in SETTINGS:
		SETTINGS['BOWTIE'] = SETTINGS['BOWTIE'][0]
	else:
		SETTINGS['BOWTIE'] = ''
	setup_bowtie()

def setup_bowtie():
	SETTINGS['cmd-bowtie-build'] = os.path.join(SETTINGS['BOWTIE'],\
			'bowtie-build')
	SETTINGS['cmd-bowtie'] = os.path.join(SETTINGS['BOWTIE'],\
			'bowtie')
	# TODO: sanity check for bowtie exists and runs!

read_config()
# change the gap markers of SeqVector globally
SeqVector.gap_markers = SETTINGS['GAP']
