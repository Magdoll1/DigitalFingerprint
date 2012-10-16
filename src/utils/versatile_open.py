import os
import sys
import gzip
import bz2

def versatile_open(file, *args):
	if file.endswith('.gz'):
		return gzip.open(file, *args)
	if file.endswith('.bz2'):
		return bz2.BZ2File(file, *args)
	return open(file, *args)

def versatile_open_func_wrap(arg_indices, kwarg_names):
	"""
	can be used to wrap a function that takes the argument "filename"
	"""
	def open_file(ref, key, cmds):
		try:
			filename = ref[key]
		except KeyError:
			return # do nothing
		if filename.endswith('.bz2'):
			os.system("bunzip2 " + filename)
			filename = filename[:-4]
			ref[key] = filename
			cmds.append("bzip2 " + filename)
		elif filename.endswith('.gz'):
			os.system("gunzip " + filename)
			filename = filename[:-3]
			ref[key] = filename
			cmds.append("gzip " + filename)

	def wrap(f):
		def g(*args, **kwargs):
			cmds = []
			args = list(args)
			for i in arg_indices: 
				open_file(args, i, cmds)
			for key in kwarg_names:
				open_file(kwargs, key, cmds)
			f(*tuple(args), **kwargs)
			for cmd in cmds:
				os.system(cmd)
		return g
	return wrap

@versatile_open_func_wrap([0],[])
def test(filename):
	print("Not doing anything, just hanging out here :D")
	with open(filename) as f:
		for line in f:
			print(line)
	
if __name__ == "__main__":
	test(sys.argv[1])
