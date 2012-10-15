import os,sys
from DigitalFingerprint.DF import DFReader
from DigitalFingerprint.SILVA import Ecoli1542_SILVA100 as ecoli


def main(df_filename):
	first_nz_pos, last_nz_pos = 1, 1542
	X = {} # df name --> ecoli coverage vec
	for df in DFReader(open(df_filename)):
			x = df.get_compressed_vec()[ecoli] # DoC of all 1,542 ecoli positions
			first_nz_pos = min(first_nz_pos, x.nonzero()[0][0]+1)
			last_nz_pos = max(last_nz_pos, x.nonzero()[0][-1]+1)
			X[df.name] = x

	f = open(df_filename + '.Ecoli_coverage.txt', 'w')
	f.write("ECOLI_POSITIONS," + ",".join(map(str, xrange(first_nz_pos, last_nz_pos+1))) + '\n')
	first_threshold_pos, last_threshold_pos = first_nz_pos, last_nz_pos
	for name, vec in X.iteritems():
		_t = vec.max() *.8
		for i in xrange(first_nz_pos-1, last_nz_pos-10):
			if vec[i:(i+10)].mean() >= _t:
				first_threshold_pos = max(first_threshold_pos, i)
				break
		for i in xrange(last_nz_pos, first_nz_pos+10, -1):
			if vec[(i-10):i].mean() >= _t:
				last_threshold_pos = min(last_threshold_pos, i)
				break
		f.write(name + ',' + ",".join(map(str, vec[(first_nz_pos-1):last_nz_pos])) + '\n')

	print >> sys.stderr, "Coverage written to {0}".format(f.name)
	print >> sys.stderr, "Use E.coli positions {0}-{1} for clustering!".format(\
			first_threshold_pos+1, last_threshold_pos+1)
	return first_threshold_pos+1, last_threshold_pos+1

if __name__ == "__main__":
	main(sys.argv[1])
		
