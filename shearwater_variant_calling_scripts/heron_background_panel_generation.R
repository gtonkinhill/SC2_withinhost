sample_list <- read.table("sample_list_1181.tsv", sep = "\t", stringsAsFactors = F, header = F)

good_samples <- readRDS("good_samples_521.rds")
normal_list_1 <- readRDS("normal_panel_list_1.rds")
normal_list_2 <- readRDS("normal_panel_list_2.rds")

run_1_pile_up <- readRDS("run_1_coverage.rds")
run_2_pile_up <- readRDS("run_2_coverage.rds")

run_1_1pct_vaf_muts <- readRDS("run_1_mutations_above_1pct_vaf.rds")
run_2_1pct_vaf_muts <- readRDS("run_2_mutations_above_1pct_vaf.rds")

shearwater_background_1 <- array(data = NA, dim = c(length(normal_list_1),29903,12))
shearwater_background_2 <- array(data = NA, dim = c(length(normal_list_2),29903,12))

background_1_mut_exc_array <- array(data = F, dim = c(29903,100))
background_2_mut_exc_array <- array(data = F, dim = c(29903,100))

#Background 1
for(i in 1:50){
  temp_counts <- run_1_pile_up[which(sample_list == normal_list_1[i]),,]
  mut_above_1pct <- run_1_1pct_vaf_muts[which(run_1_1pct_vaf_muts$sampleID == normal_list_1[i]),]

  background_1_mut_exc_array[as.numeric(mut_above_1pct$pos),i] <- T
  
  temp_counts[background_1_mut_exc_array[,i],] <- 0
  
  for(k in 1:12){
    shearwater_background_1[i,,k] <- temp_counts[,k]
  }
}

for(i in 51:100){
  temp_counts <- run_2_pile_up[which(sample_list == normal_list_1[i]),,]
  mut_above_1pct <- run_2_1pct_vaf_muts[which(run_2_1pct_vaf_muts$sampleID == normal_list_1[i]),]
  
  background_1_mut_exc_array[as.numeric(mut_above_1pct$pos),i] <- T
  
  temp_counts[background_1_mut_exc_array[,i],] <- 0
  
  for(k in 1:12){
    shearwater_background_1[i,,k] <- temp_counts[,k]
  }
}
saveRDS(shearwater_background_1,"background_panel_1.rds")

#Background_2
for(i in 1:50){
  temp_counts <- run_1_pile_up[which(sample_list == normal_list_2[i]),,]
  mut_above_1pct <- run_1_1pct_vaf_muts[which(run_1_1pct_vaf_muts$sampleID == normal_list_2[i]),]
  
  background_2_mut_exc_array[as.numeric(mut_above_1pct$pos),i] <- T
  
  temp_counts[background_2_mut_exc_array[,i],] <- 0
  
  for(k in 1:12){
    shearwater_background_2[i,,k] <- temp_counts[,k]
  }
}

for(i in 51:100){
  temp_counts <- run_2_pile_up[which(sample_list == normal_list_2[i]),,]
  mut_above_1pct <- run_2_1pct_vaf_muts[which(run_2_1pct_vaf_muts$sampleID == normal_list_2[i]),]
  
  background_2_mut_exc_array[as.numeric(mut_above_1pct$pos),i] <- T
  
  temp_counts[background_2_mut_exc_array[,i],] <- 0
  
  for(k in 1:12){
    shearwater_background_2[i,,k] <- temp_counts[,k]
  }
}
saveRDS(shearwater_background_2,"background_panel_2.rds")

background_exc_count <- as.data.frame(array(data = NA, dim = c(29903,3)))
colnames(background_exc_count) <- c("pos","exc_bg_1","exc_bg_2")
background_exc_count$pos <- 1:29903
for(i in 1:29903){
  background_exc_count[i,"exc_bg_1"] <- length(which(background_1_mut_exc_array[i,]))
  background_exc_count[i,"exc_bg_2"] <- length(which(background_2_mut_exc_array[i,]))
}

ggplot(background_exc_count) +
  geom_point(aes(x = pos, y = exc_bg_1), size = 0.5) +
  theme_bw() +
  labs(x = "Position", y = "Panel 1 exc. due to VAF above 1%")
ggsave("background_panel_1_exc_sites.pdf", width = 15, height = 5)

ggplot(background_exc_count) +
  geom_point(aes(x = pos, y = exc_bg_2), size = 0.5) +
  theme_bw() +
  labs(x = "Position", y = "Panel 2 exc. due to VAF above 1%")
ggsave("background_panel_2_exc_sites.pdf", width = 15, height = 5)

ggplot(background_exc_count) +
  geom_point(aes(x = pos, y = exc_bg_1), size = 0.5, color = "red", alpha = 0.5) +
  geom_point(aes(x = pos, y = exc_bg_2), size = 0.5, color = "blue", alpha = 0.5) +
  theme_bw() +
  labs(x = "Position", y = "Excluded due to VAF above 1%")
ggsave("background_panel_1_and_2_exc_sites.pdf", width = 15, height = 5)

