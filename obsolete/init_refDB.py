from settings import *

def process_gapped_fasta(filename, output_dir):
	"""
	Given a gapped fasta file, in output_dir,
	make a copy (in DNA), and a .ungapped version,
	and a .gap_map mapping file
	"""
	start_t = datetime.datetime.now()
	_basename = os.path.basename(filename)
	f_copy = open(os.path.join(output_dir, _basename), 'w')
	f_ungapped = open(os.path.join(output_dir, _basename+'.ungapped'), 'w')
	f_gapmap = open(os.path.join(output_dir, _basename+'.gap_map'), 'w')

	msg_len = 0
	with open(filename) as f:
		for r in SeqIO.parse(f, 'fasta'):
			msg = "Processing seq "  + r.id
			sys.stdout.write("\b"*msg_len + msg) 
			msg_len = len(msg)

			seq = str(r.seq).replace('U', 'T')
			f_copy.write(">{id}\n{seq}\n".format(id=r.id, seq=seq))

			ungapped_pos = []
			ungapped_seq = ''
			for i,s in enumerate(seq):
				if s not in SETTINGS['GAP']:
					ungapped_seq += s
					ungapped_pos.append(str(i))

			f_ungapped.write(">{id}\n{seq}\n".format(\
					id=r.id, seq=ungapped_seq))
			f_gapmap.write("{id}\t{map}\n".format(\
					id=r.id, map=",".join(ungapped_pos)))
	end_t = datetime.datetime.now()

	f_copy.close()
	f_ungapped.close()
	f_gapmap.close()
	with open(SETTINGS['LOG'], 'a+') as f_log:
		f_log.write(time.ctime() + '\n')
		f_log.write("initDB: processed gapped fasta file {0}\n".format(filename))
		f_log.write("Output stored in {0}\n".format(output_dir))
		f_log.write("Ungapped version: {0}\n".format(f_ungapped.name))
		f_log.write("Gap map file: {0}\n".format(f_gapmap.name))
		f_log.write("Processing time: {0}\n".format(end_t-start_t))

if __name__ == "__main__":
	parser = OptionParser()
	parser.add_option("-d", "--output-dir", dest="output_dir",\
			default=SETTINGS['REF'], \
			help="output directory (default set by config file)")
	parser.add_option("-i", "--input-fasta", dest="input_fasta",\
			help="input (gapped fasta) filename")

	options, args = parser.parse_args()
	process_gapped_fasta(options.input_fasta, options.output_dir)

