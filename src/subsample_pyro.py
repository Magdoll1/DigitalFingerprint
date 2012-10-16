import os, sys, glob
from Pyro import Pyro, PyroSampler
from DF import DFWriter

if __name__ == "__main__":
	pattern = sys.argv[1]
	subsample_size = int(sys.argv[2])
	output = sys.argv[3]

	f = open(output, 'w')
	w = DFWriter(f)
	for file in glob.iglob(pattern):
		pyro = Pyro(os.path.basename(file), file, None, False)
		sampler = PyroSampler(pyro)
		df = sampler.subsample(subsample_size)
		w.write(df)
	f.close()
