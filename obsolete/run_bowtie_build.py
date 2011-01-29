from settings import *

def run_bowtie_build(filename):
	seed = int(time.time())
	out_log = filename + '.bowtie-build.log'
	cmd = "{cmd} -f --seed {seed} {input} {out} > {out_log}".format(\
			cmd=SETTINGS['cmd-bowtie-build'],\
			seed=seed,\
			input=filename,\
			out=filename,\
			out_log=out_log)
	print >> sys.stderr, "running command", cmd
	start_t = datetime.datetime.now()
	os.system(cmd)
	end_t = datetime.datetime.now()
	with open(SETTINGS['LOG'], 'a+') as f_log:
		f_log.write(time.ctime() + '\n')
		f_log.write("Ran bowtie-build for {0}\n".format(filename))
		f_log.write("cmd: {0}\n".format(cmd))
		f_log.write("seed: {0}\n".format(seed))
		f_log.write("bowtie-build running log: {0}\n".format(out_log))
		f_log.write("Processing time: {0}\n".format(end_t-start_t))

if __name__ == "__main__":
	parser = OptionParser()
	parser.add_option("-i", "--input-fasta", dest="input_fasta")

	options, args = parser.parse_args()
	run_bowtie_build(options.input_fasta)

