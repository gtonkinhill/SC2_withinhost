args = commandArgs(TRUE)
run_number = args[1]
sample_name = args[2]
bam_file = args[3]

library("GenomicRanges")
library("deepSNV")
library("Rsamtools")

nucleotides = c("A", "T", "C", "G", "-", "INS", "a", "t", "c", "g", "_", "ins")

minPhred=30
samFlag=3844
mapq=30

test.matrix <- matrix(0, ncol = length(nucleotides), nrow = 29903)
colnames(test.matrix) <- nucleotides

test.matrix = bam2R(bam_file, "MN908947.3", 1, 29903, q=minPhred, mask=samFlag, mq=mapq)[, nucleotides]

mode(test.matrix) = "integer"

write.table(test.matrix, file = paste0("run_",run_number,"/",sample_name,"_run",run_number,"_count.csv"), sep = ",",quote = F,row.names = F,col.names = T)
