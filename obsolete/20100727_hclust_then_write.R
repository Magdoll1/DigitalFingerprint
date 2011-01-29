for (e in commandArgs()) {
	ta <- strsplit(e, '=');
	if (ta[[1]][1] == '--input') {
		input_filename <- ta[[1]][2]; 
	}
	if (ta[[1]][1] == '--output') {
		output_filename <- ta[[1]][2];
	}
}

x <- read.table(input_filename, sep=',', row.names=1);
c <- hclust( dist( x ) );
write.table( c$merge, output_filename, append=TRUE );

