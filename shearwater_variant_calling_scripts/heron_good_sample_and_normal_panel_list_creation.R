sample_info <-  read.table("sample_list_1181.tsv", sep = "\t", stringsAsFactors = F, header = F)
colnames(sample_info) <- "sampleID"

im3_calculation <- read.table("initial_analyses_using_im3_mut_read_above_500x_table/Heron_all_mismatches.txt", sep = "\t", stringsAsFactors = F, header = T)
replicate_info <- read.table("replicate_meta.tsv", sep = "\t", stringsAsFactors = F, header = T)

for(i in 1:nrow(sample_info)){
  sample_info$Ct[i] <- unique(im3_calculation$Ct[which(im3_calculation$sampleid == sample_info$sampleID[i])])
  sample_info$rho[i] <- unique(im3_calculation$rho[which(im3_calculation$sampleid == sample_info$sampleID[i])])
  sample_info$total_samples[i] <- replicate_info$n_samples[which(replicate_info$sample_id == sample_info$sampleID[i])]
  sample_info$donor_id[i] <- replicate_info$biosample_sources[which(replicate_info$sample_id == sample_info$sampleID[i])] 
}

sample_info$seq_samples <- NA

for(i in 1:nrow(sample_info)){
  if(sample_info$donor_id[i] != ""){
    sample_info$seq_samples[i] <- length(which(sample_info$donor_id == sample_info$donor_id[i]))
  }
}

sample_coverage_summary <- readRDS("summary_coverage_statistics_per_sample.rds")

for(i in 1:nrow(sample_info)){
  sample_info$run_1_prop_above_500x[i] <- sample_coverage_summary$run_1_prop_above_500x[which(sample_coverage_summary$sampleID == sample_info$sampleID[i])]
  sample_info$run_2_prop_above_500x[i] <- sample_coverage_summary$run_2_prop_above_500x[which(sample_coverage_summary$sampleID == sample_info$sampleID[i])]
  sample_info$run_1_median_coverage[i] <- sample_coverage_summary$run_1_50pct[which(sample_coverage_summary$sampleID == sample_info$sampleID[i])]
  sample_info$run_2_median_coverage[i] <- sample_coverage_summary$run_2_50pct[which(sample_coverage_summary$sampleID == sample_info$sampleID[i])]
}

ggplot(sample_info[which(sample_info$run_2_prop_above_500x > 0.25),]) +
  geom_point(aes(x = run_1_prop_above_500x, y = rho)) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  labs(x = "Run 1 - Proportion of genome with coverage above 500x", y = "Rho")
ggsave("run_1_genome_prop_above_500x_vs_rho.pdf", width = 10, height = 10)

ggplot(sample_info) +
  geom_point(aes(x = run_2_prop_above_500x, y = rho)) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  labs(x = "Run 2 - Proportion of genome with coverage above 500x", y = "Rho")
ggsave("run_2_genome_prop_above_500x_vs_rho.pdf", width = 10, height = 10)

bad_samples_fail_rho <- sample_info$sampleID[which(sample_info$rho > 0.02)]
saveRDS(bad_samples_fail_rho, "bad_samples_fail_rho_492.rds")
bad_samples_fail_90pct_500x_dep <- sample_info$sampleID[which(sample_info$run_1_prop_above_500x < 0.9 | sample_info$run_2_prop_above_500x < 0.9)]
saveRDS(bad_samples_fail_90pct_500x_dep, "bad_samples_fail_90pct_above_500x_both_runs_588.rds")

good_samples <- sample_info$sampleID[which((sample_info$rho <= 0.02) & (sample_info$run_1_prop_above_500x >= 0.9 & sample_info$run_2_prop_above_500x >= 0.9) & (abs(sample_info$run_1_median_coverage - sample_info$run_2_median_coverage) <= 20000))]
good_samples_non_dup <- sample_info$sampleID[which((sample_info$rho <= 0.02) & (sample_info$run_1_prop_above_500x >= 0.9 & sample_info$run_2_prop_above_500x >= 0.9) & (abs(sample_info$run_1_median_coverage - sample_info$run_2_median_coverage) <= 20000) & (is.na(sample_info$seq_samples) | (sample_info$seq_samples == 1)))]

combined_normal_panel <- sample(good_samples_non_dup, 200, replace = F)
normal_panel_1 <- combined_normal_panel[1:100]
normal_panel_2 <- combined_normal_panel[101:200]

saveRDS(good_samples, "good_samples_521.rds")
saveRDS(normal_panel_1, "normal_panel_list_1.rds")
saveRDS(normal_panel_2, "normal_panel_list_2.rds")