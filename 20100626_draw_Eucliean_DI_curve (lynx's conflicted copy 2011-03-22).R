for (e in commandArgs()) {
    ta <- strsplit(e, '=');
    if (ta[[1]][1] == '--input') {
        input_filename <- ta[[1]][2];
    }
    if (ta[[1]][1] == '--output') {
        output_filename <- ta[[1]][2];
    }
}

filename <- input_filename;
png(output_filename);
f <- file(filename);
sample_name <- readLines(f, n=1)[1];
r <- read.table(filename, row.names=1, sep=',', skip=1);
y <- c();
num <- length(row.names(r));
ymax <- 0;
for (i in 1:num) { y <- c(y, sum(r[i,])/length(r[i,])); ymax <- max(ymax, max(r[i,])); }
ymax <- (round(ymax / 0.2)+1) * 0.2;
par(mar=par()$mar+c(2,0,0,0));
plot(y,ylim=c(0,ymax),axes=FALSE,xlab='',ylab='Distance to full DF');
for (i in 1:num) { arrows(i, min(r[i,]), i, max(r[i,]), length=0); };
axis(1,1:length(y),row.names(r),las=2);
axis(2,seq(0,ymax,0.2),las=1);
# below two lines for drawing blue horizontal lines for visual aid
#library(fields);
#yline(seq(0,ymax,0.2),col=4);
title(sample_name);
mtext("Subsample size",side=1,line=5);
dev.off()
