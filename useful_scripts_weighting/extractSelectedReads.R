library(data.table)

INPUTFILE="output.paf" #mapping output in paf format
READIDS="readids.repetitive.txt" #read id, one per line
OUTPUTFILE="output.repetitive.paf" #filtered paf output with specified read ids

allMappings = read.table(INPUTFILE, header=F, fill=TRUE)
selectedReadIds = read.table(READIDS, header=F, fill=TRUE)
selectedMappings = allMappings[allMappings$V1 %in% selectedReadIds$V1,] 
write.table(selectedMappings, OUTPUTFILE, quote = FALSE, sep = "\t", row.names = F, col.names = F)

# We could have used "grep -f" for the above task, but grep is terribly slow for this op.
