import os, sys, glob
from DigitalFingerprint.DF import DFWriter
from DigitalFingerprint.Pyro import Pyro

if __name__ == "__main__":
	from optparse import OptionParser

	parser = OptionParser()
	parser.add_option("-p", "--pattern", dest="pattern", help="pattern for list of input files")
	parser.add_option("-o", "--output", dest="output", help="output filename (suffix should be .DF)")
	parser.add_option("--dup-ok", dest="rename_dups", action="store_true", default=False, help="ok to have duplicate IDs (renames them)")

	options, args = parser.parse_args()
	pattern, output = options.pattern, options.output

	f = open(output, 'w')
	w = DFWriter(f)
	for file in glob.iglob(pattern):
		print >> sys.stderr, "processing", file
		name = os.path.basename(file)
		pyro = Pyro(name, file, None, options.rename_dups)
		df = pyro.make_DF()
		w.write(df)
	f.close()

