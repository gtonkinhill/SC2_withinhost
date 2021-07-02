run_1_coverage <- readRDS("run_1_coverage.rds")
run_2_coverage <- readRDS("run_2_coverage.rds")

sample_coverage_summary <- as.data.frame(array(data = NA, dim = c(1181,19)))
colnames(sample_coverage_summary) <- c("sampleID",
                                       "run_1_mean","run_1_min","run_1_max","run_1_10pct","run_1_25pct","run_1_50pct","run_1_75pct","run_1_90pct","run_1_prop_above_500x",
                                       "run_2_mean","run_2_min","run_2_max","run_2_10pct","run_2_25pct","run_2_50pct","run_2_75pct","run_2_90pct","run_2_prop_above_500x")

sample_list <- read.table("sample_list_1181.tsv", sep = "\t", stringsAsFactors = F, header = F)

sample_coverage_summary$sampleID <- sample_list

run_1_normalised_coverage <- array(data = NA, dim = c(29903,1181))
run_2_normalised_coverage <- array(data = NA, dim = c(29903,1181))

for(i in 1:1181){
  sample_coverage_summary$run_1_mean[i] <- sum(run_1_coverage[i,,]) / 29903
  sample_coverage_summary$run_1_min[i] <- min(rowSums(run_1_coverage[i,,]))
  sample_coverage_summary$run_1_max[i] <- max(rowSums(run_1_coverage[i,,]))
  sample_coverage_summary$run_1_10pct[i] <- quantile(rowSums(run_1_coverage[i,,]), 0.1)
  sample_coverage_summary$run_1_25pct[i] <- quantile(rowSums(run_1_coverage[i,,]), 0.25)
  sample_coverage_summary$run_1_50pct[i] <- quantile(rowSums(run_1_coverage[i,,]), 0.5)
  sample_coverage_summary$run_1_75pct[i] <- quantile(rowSums(run_1_coverage[i,,]), 0.75)
  sample_coverage_summary$run_1_90pct[i] <- quantile(rowSums(run_1_coverage[i,,]), 0.9)
  sample_coverage_summary$run_1_prop_above_500x[i] <- length(which(rowSums(run_1_coverage[i,,]) >= 500)) / 29903
  
  run_1_normalised_coverage[,i] <- rowSums(run_1_coverage[i,,]) / sample_coverage_summary$run_1_mean[i]
  
  sample_coverage_summary$run_2_mean[i] <- sum(run_2_coverage[i,,]) / 29903
  sample_coverage_summary$run_2_min[i] <- min(rowSums(run_2_coverage[i,,]))
  sample_coverage_summary$run_2_max[i] <- max(rowSums(run_2_coverage[i,,]))
  sample_coverage_summary$run_2_10pct[i] <- quantile(rowSums(run_2_coverage[i,,]), 0.1)
  sample_coverage_summary$run_2_25pct[i] <- quantile(rowSums(run_2_coverage[i,,]), 0.25)
  sample_coverage_summary$run_2_50pct[i] <- quantile(rowSums(run_2_coverage[i,,]), 0.5)
  sample_coverage_summary$run_2_75pct[i] <- quantile(rowSums(run_2_coverage[i,,]), 0.75)
  sample_coverage_summary$run_2_90pct[i] <- quantile(rowSums(run_2_coverage[i,,]), 0.9)
  sample_coverage_summary$run_2_prop_above_500x[i] <- length(which(rowSums(run_2_coverage[i,,]) >= 500)) / 29903
  
  run_2_normalised_coverage[,i] <- rowSums(run_2_coverage[i,,]) / sample_coverage_summary$run_2_mean[i]
}

colnames(run_1_normalised_coverage) <- sample_list$V1
colnames(run_2_normalised_coverage) <- sample_list$V1

saveRDS(object = sample_coverage_summary, file = "summary_coverage_statistics_per_sample.rds")
saveRDS(object = run_1_normalised_coverage, file = "run_1_normalised_coverage_per_sample.rds")
saveRDS(object = run_2_normalised_coverage, file = "run_2_normalised_coverage_per_sample.rds")

saveRDS(object = sample_coverage_summary, file = "summary_coverage_statistics_per_sample.rds")
saveRDS(object = run_1_normalised_coverage, file = "run_1_normalised_coverage_per_sample.rds")
saveRDS(object = run_2_normalised_coverage, file = "run_2_normalised_coverage_per_sample.rds")

library(tidyverse)

sample_coverage_summary$run_1_med_idx[order(sample_coverage_summary$run_1_50pct)] <- 1:1181
sample_coverage_summary$run_2_med_idx[order(sample_coverage_summary$run_2_50pct)] <- 1:1181
sample_coverage_summary$run_1_prop_above_500_idx[order(sample_coverage_summary$run_1_prop_above_500x)] <- 1:1181

ggplot(sample_coverage_summary) +
  geom_point(aes(x = run_1_med_idx, y = run_1_50pct, color = "Run 1 - Median")) +
  geom_point(aes(x = run_1_med_idx, y = run_1_25pct, color = "Run 1 - Lower Quartile")) +
  geom_point(aes(x = run_1_med_idx, y = run_1_75pct, color = "Run 1 - Upper Quartile")) +
  theme_bw() +
  labs(x = "Sample number (ordered by Run 1 - Median)", y = "Coverage", color = "Legend")
ggsave(filename = "run_1_median_quartile_coverage.pdf", width = 10, height = 5)

ggplot(sample_coverage_summary) +
  geom_point(aes(x = run_2_med_idx, y = run_2_50pct, color = "Run 2 - Median")) +
  geom_point(aes(x = run_2_med_idx, y = run_2_25pct, color = "Run 2 - Lower Quartile")) +
  geom_point(aes(x = run_2_med_idx, y = run_2_75pct, color = "Run 2 - Upper Quartile")) +
  theme_bw() +
  labs(x = "Sample number (ordered by Run 2 - Median)", y = "Coverage", color = "Legend")
ggsave(filename = "run_2_median_quartile_coverage.pdf", width = 10, height = 5)

ggplot(sample_coverage_summary) +
  geom_point(aes(x = run_1_50pct, y = run_2_50pct)) +
  geom_abline(slope = 1, intercept = 0, color = "red", lty = 2) +
  theme_bw() +
  labs(x = "Run 1 - Median Coverage", y = "Run 2 - Median Coverage") +
  scale_y_continuous(limits = c(0,50000)) +
  scale_x_continuous(limits = c(0,150000))
ggsave(filename = "run_1_vs_run_2_median_coverage.pdf", width = 15, height = 5)

ggplot(sample_coverage_summary) +
  geom_point(aes(x = run_1_prop_above_500x, y = run_2_prop_above_500x)) +
  scale_y_continuous(limits = c(0,1)) +
  scale_x_continuous(limits = c(0,1)) +
  theme_bw() +
  labs(x = "Run 1 - Proportion of genome above 500x coverage", y = "Run 2 - Proportion of genome above 500x coverage")
ggsave(filename = "run_1_vs_run_2_genome_prop_above_500x_coverage.pdf", width = 10, height = 10)

per_site_norm_coverage <- as.data.frame(array(data = NA, dim = c(29903,17)))
colnames(per_site_norm_coverage) <- c("position",
                                       "run_1_mean","run_1_min","run_1_max","run_1_10pct","run_1_25pct","run_1_50pct","run_1_75pct","run_1_90pct",
                                       "run_2_mean","run_2_min","run_2_max","run_2_10pct","run_2_25pct","run_2_50pct","run_2_75pct","run_2_90pct")

per_site_norm_coverage$position <- 1:29903
for(i in 1:29903){
  per_site_norm_coverage$run_1_mean[i] <- mean(run_1_normalised_coverage[i,])
  per_site_norm_coverage$run_1_min[i] <- min(run_1_normalised_coverage[i,])
  per_site_norm_coverage$run_1_max[i] <- max(run_1_normalised_coverage[i,])
  per_site_norm_coverage$run_1_10pct[i] <- quantile(run_1_normalised_coverage[i,], 0.1)
  per_site_norm_coverage$run_1_25pct[i] <- quantile(run_1_normalised_coverage[i,], 0.25)
  per_site_norm_coverage$run_1_50pct[i] <- quantile(run_1_normalised_coverage[i,], 0.5)
  per_site_norm_coverage$run_1_75pct[i] <- quantile(run_1_normalised_coverage[i,], 0.75)
  per_site_norm_coverage$run_1_90pct[i] <- quantile(run_1_normalised_coverage[i,], 0.9)
  
  per_site_norm_coverage$run_2_mean[i] <- mean(run_2_normalised_coverage[i,], na.rm = T)
  per_site_norm_coverage$run_2_min[i] <- min(run_2_normalised_coverage[i,], na.rm = T)
  per_site_norm_coverage$run_2_max[i] <- max(run_2_normalised_coverage[i,], na.rm = T)
  per_site_norm_coverage$run_2_10pct[i] <- quantile(run_2_normalised_coverage[i,], 0.1, na.rm = T)
  per_site_norm_coverage$run_2_25pct[i] <- quantile(run_2_normalised_coverage[i,], 0.25, na.rm = T)
  per_site_norm_coverage$run_2_50pct[i] <- quantile(run_2_normalised_coverage[i,], 0.5, na.rm = T)
  per_site_norm_coverage$run_2_75pct[i] <- quantile(run_2_normalised_coverage[i,], 0.75, na.rm = T)
  per_site_norm_coverage$run_2_90pct[i] <- quantile(run_2_normalised_coverage[i,], 0.9, na.rm = T)
}

heron_primers <- read.table("nCoV-2019.bed", sep = "\t", stringsAsFactors = F, header = F)
primer_positions <- as.data.frame(array(data = NA, dim = c(98,3)))
colnames(primer_positions) <- c("primerID","start","end")

for(i in 1:nrow(primer_positions)){
  if(!(i %in% c(7,9,14,15,18,21,44,45,46,76,89))){
    primer_positions$primerID[i] <- paste0("nCoV-2019_",i)
    primer_positions$start[i] <- heron_primers$V3[which(heron_primers$V4 == paste0("nCoV-2019_",i,"_LEFT"))] + 1
    primer_positions$end[i] <- heron_primers$V2[which(heron_primers$V4 == paste0("nCoV-2019_",i,"_RIGHT"))]
  } else{
    left_name <- heron_primers$V4[grepl(paste0("nCoV-2019_",i,"_LEFT_alt"), x = heron_primers$V4)]
    right_name <- heron_primers$V4[grepl(paste0("nCoV-2019_",i,"_RIGHT_alt"), x = heron_primers$V4)]
    primer_positions$primerID[i] <- paste0("nCoV-2019_",i,"_alt")
    primer_positions$start[i] <- heron_primers$V3[which(heron_primers$V4 == left_name)] + 1
    primer_positions$end[i] <- heron_primers$V3[which(heron_primers$V4 == right_name)]
  }
}

per_site_norm_coverage$primers_covered <- ""
for(j in 1:nrow(primer_positions)){
  for(k in primer_positions$start[j]:primer_positions$end[j]){
    per_site_norm_coverage$primers_covered[k] <- paste0(per_site_norm_coverage$primers_covered[k],primer_positions$primerID[j],";") 
  }
}

for(i in 1:nrow(per_site_norm_coverage)){
  per_site_norm_coverage$num_primers[i] <- str_count(per_site_norm_coverage$primers_covered[i],";") 
}

ggplot(per_site_norm_coverage) +
  geom_point(aes(x = position, y = run_1_50pct, color = as.character(num_primers)), size = 0.2) +
  theme_bw() +
  labs(x = "Genomic position", y = "Median of normalised coverage across Run 1 samples", color = "Number of amplicons")
ggsave("run_1_genomic_pos_median_normalised_coverage_num_amplicons.pdf", height = 5, width = 15)

ggplot(per_site_norm_coverage) +
  geom_point(aes(x = position, y = run_2_50pct), size = 0.2) +
  theme_bw() +
  labs(x = "Genomic position", y = "Median of normalised coverage across Run 2 samples")
ggsave("run_2_genomic_pos_median_normalised_coverage.pdf", height = 5, width = 15)

ggplot(per_site_norm_coverage) +
  geom_point(aes(x = position, y = run_1_50pct, color = "Run 1"), size = 0.2) +
  geom_point(aes(x = position, y = run_2_50pct, color = "Run 2"), size = 0.2) +
  theme_bw() +
  labs(x = "Genomic position", y = "Median of normalised coverage across samples in a run", color = "Legend")
ggsave("both_runs_genomic_pos_median_normalised_coverage.pdf", height = 5, width = 15)
  
ggplot(per_site_norm_coverage) +
  geom_point(data = per_site_norm_coverage[c(1:385,665:29903),], aes(x = run_1_50pct, y = run_2_50pct, color = "Outside amplicon 2")) +
  geom_point(data = per_site_norm_coverage[386:664,], aes(x = run_1_50pct, y = run_2_50pct, color = "In amplicon 2")) +
  theme_bw() +
  scale_x_continuous(limits = c(0,4.5)) +
  scale_y_continuous(limits = c(0,4.5)) +
  labs(x = "Median of normalised coverage at each position across samples in Run 1", y = "Median of normalised coverage at each position across samples in Run 2", color = "Legend")
ggsave("run_1_vs_run_2_genome_pos_median_normalised_coverage.pdf", width = 10, height = 10)
