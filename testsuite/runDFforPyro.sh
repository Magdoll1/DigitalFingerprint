#!/bin/bash

REF_FILE=silva.bacteria/silva.bacteria.fasta
script_DIR=../scripts/
INPUT_DIR=16Spyro_28F_200seq

if [ ! -f $REF_FILE ]; then
	echo "Downloading SILVA reference seed alignment from Mothur website...."
	wget http://www.mothur.org/w/images/9/98/Silva.bacteria.zip
	unzip Silva.bacteria.zip
	REF_FILE=silva.bacteria/silva.bacteria.fasta
fi

if [ -d $INPUT_DIR ]; then
	base=`basename $INPUT_DIR`
#	for FILE in `ls $base/*.fna`
##	do 
#		mothur "#align.seqs(candidate=$FILE, template=$REF_FILE, processors=4, flip=F)"
#	done
	python $script_DIR/Pyro_process.py -p "$base/*.align" --dup-ok -o $base/$base.DF
	echo "DF file written to $base/$base.DF"
	echo "Already know recommended E.coli range is 20-342. Running clustering...."
	python $script_DIR/run_clustering.py -f $base/$base.DF -r 20,342 -i Simpson -d $base/$base.Simpson20to342.DI.txt \
		-o $base/$base.Simpson20to342.tree
else
	echo "$INPUT_DIR is not a directory! Terminate."
fi

