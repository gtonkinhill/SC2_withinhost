sample_list <- read.table("sample_list_1181.tsv", sep = "\t", stringsAsFactors = F, header = F)
bam_files <- vector(length = nrow(sample_list))
run_number <- 2

output_files <- list.files(paste0("Data/run_",run_number))

for(i in 1:nrow(sample_list)){
  bam_path <- paste0("replicates_experiment/input/run_",run_number,"/",sample_list[i,])
  bam_files[i] <-  paste0(bam_path,"/",list.files(bam_path)[grepl(pattern = "primertrimmed.sorted.bam$", list.files(bam_path))])

  if(!(paste0("Data/run_",run_number,"/",sample_list[i,],"_run",run_number,"_count.csv") %in% output_files)){
    system(sprintf("bsub -q normal -M5000 -R'select[mem>5000] rusage[mem=5000] span[hosts=1]' -n1 -e Data/run_%s/logs/%s.err -o Data/run_%s/logs/%s.out /software/R-3.6.1/bin/Rscript al28_code/heron_genome_wide_pileup.R %s %s %s",
                   run_number, sample_list[i,], run_number, sample_list[i,], run_number, sample_list[i,], bam_files[i]))
  }
}