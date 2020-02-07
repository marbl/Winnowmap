#Purpose: Draw a density plot indicating minimizer density across a chromosome
#Apply hard masking cutoff and see how it affects minimizer distribution
#Assumes a file called "minimizers.txt" containing minimizer locations on a contguous sequence

library(data.table)
library(plyr)

args <- commandArgs(trailingOnly = TRUE)

#MAX_OCC should be +ve positive value to enable filtering
#OR set it to -1 to disable filtering
MAX_OCC=as.numeric(args[1])	

#Reference sequence id assigned by minimap2 to that chromosome
#this should be the serial no. of that sequence (0-based)
#you can check .fai file 
REFSEQID=as.numeric(args[2])
INPUTFILE="minimizers.txt"
OUTPUTFILE="density.pdf" 

allMinimizers = fread(INPUTFILE, header=F)

if (MAX_OCC >= 0)
{
	counts = count(allMinimizers, 'V3')
	filteredMinimizers = subset(counts, counts$freq <= MAX_OCC)
	preservedMinimizers = allMinimizers[allMinimizers$V3 %in% filteredMinimizers$V3]
	chrXMinimizers = subset(preservedMinimizers, preservedMinimizers$V1 == REFSEQID)
} else { 
	chrXMinimizers = subset(allMinimizers, allMinimizers$V1 == REFSEQID)
}

pdf(OUTPUTFILE, width=6,height=4)
d <- density(chrXMinimizers$V2, bw = 500)
plot(d, main="Minimizer density across chromosome X")
dev.off()
