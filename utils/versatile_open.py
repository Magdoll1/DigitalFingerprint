import gzip
import bz2

def versatile_open(file, *args):
	if file.endswith('.gz'):
		return gzip.open(file, *args)
	if file.endswith('.bz2'):
		return bz2.BZ2File(file, *args)
	return open(file, *args)

