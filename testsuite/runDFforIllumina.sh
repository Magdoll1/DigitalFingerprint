#!/bin/bash

DF_DIR=../
INPUT_DIR=16Sillumina_QinTest

bunzip2 $INPUT_DIR/*.bowtied.bz2

if [ -d $INPUT_DIR ]; then
	base=`basename $INPUT_DIR`
#	python $DF_DIR/Read_processbowtie.py -p "$base/*.bowtied" -o $base/$base.DF
	echo "DF file written to $base/$base.DF"
	echo "Already know recommended E.coli range is 641-1356. Running clustering...."
	python $DF_DIR/clustering.py -f $base/$base.DF -r 641,1356 -i Entropy \
		-d $base/$base.Entropy641to1356.DI.txt -o $base/$base.Entropy641to1456.tree
else
	echo "$INPUT_DIR is not a directory! Terminate."
fi

