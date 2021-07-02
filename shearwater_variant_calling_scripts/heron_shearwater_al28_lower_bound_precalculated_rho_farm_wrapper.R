sample_list <- read.table("sample_list_1181.tsv", sep = "\t", stringsAsFactors = F, header = F)

background_panels <- c(1,2)
run_list <- c(1,2)
sample_start <- seq(from = 1, to = 1151, by = 50)
sample_end <- c(seq(from = 50, to = 1150, by = 50),1181)

shearwater_fun <- "heron_shearwater_al28_lower_bound_precalculated_rho_farm_script.R"

for(i in 1:length(background_panels)){
  for(j in 1:length(run_list)){
    for(k in 1:length(sample_start)){
      job_id <- paste0("b",background_panels[i],"_r",run_list[j],"_s",sample_start[k],"_",sample_end[k])
      log_path <- paste0("logs/",job_id)
      system(sprintf("bsub -J %s -q normal -R 'select[mem>=5000] rusage[mem=5000]' -M5000 -e %s.err -o %s.out /software/R-3.6.1/bin/Rscript %s %s %s %s %s", job_id, log_path, log_path, shearwater_fun, background_panels[i], run_list[j], sample_start[k], sample_end[k]))
    }
  }
}