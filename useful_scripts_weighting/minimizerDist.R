#Purpose: Draw a histogram indicating minimizer density across a chromosome
#Apply hard masking cutoff and see how it affects minimizer distribution

library(data.table)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)

MAX_OCC=as.numeric(args[1])
INPUTFILE="minimizers.txt"
OUTPUTFILE="hist.pdf" 

allMinimizers = fread(INPUTFILE, header=F)
counts = count(allMinimizers, 'V2')
filteredMinimizers = subset(counts, counts$freq <= MAX_OCC)
preservedMinimizers = allMinimizers[allMinimizers$V2 %in% filteredMinimizers$V2]

pdf(OUTPUTFILE, width=6,height=4)
hist(preservedMinimizers$V1, breaks = 300, xlab="chromosome X", ylab="count of minimizers", main="Minimizer density across chromosome X")
dev.off()
