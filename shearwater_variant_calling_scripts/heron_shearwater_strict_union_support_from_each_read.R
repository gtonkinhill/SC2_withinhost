sw_1_1 <- readRDS("shearwater_muts_all_samples_qval0.05_mincov100x_deltaVAF0.05merge_panel_1_run_1.rds")
sw_1_2 <- readRDS("shearwater_muts_all_samples_qval0.05_mincov100x_deltaVAF0.05merge_panel_1_run_2.rds")

initial_sw_1_1 <- readRDS("initial_shearwater_calls_background_panel_1_run_1.rds")
initial_sw_1_2 <- readRDS("initial_shearwater_calls_background_panel_1_run_2.rds")

run_1_coverage <- readRDS("run_1_coverage.rds")
run_2_coverage <- readRDS("run_2_coverage.rds")

ref_list <- c("A","T","C","G","-","INS")

sample_list <- read.table("sample_list_1181.tsv",sep = "\t")
sample_list <- as.vector(t(sample_list))

initial_sw_1_1$mut_site <- paste(initial_sw_1_1$pos,initial_sw_1_1$ref,initial_sw_1_1$mut,sep = "_")
initial_sw_1_1$mut_id <- paste(initial_sw_1_1$sampleID,initial_sw_1_1$mut_site,sep = "_")

initial_sw_1_2$mut_site <- paste(initial_sw_1_2$pos,initial_sw_1_2$ref,initial_sw_1_2$mut,sep = "_")
initial_sw_1_2$mut_id <- paste(initial_sw_1_2$sampleID,initial_sw_1_2$mut_site, sep = "_")

call_list <- readRDS("shearwater_intermediate_intersect_mut_list_deltaVAF0.05merge.rds")

calls <- as.data.frame(array(data = NA, dim = c(length(call_list),21)))
colnames(calls) <- c("mut_id","sample","pos","ref","mut",
                           "r1_xfw","r1_xbw","r1_nfw","r1_nbw","r1_vaf","r1_pval","r1_strict","r1_qval",
                           "r2_xfw","r2_xbw","r2_nfw","r2_nbw","r2_vaf","r2_pval","r2_strict","r2_qval")
for(i in 1:length(call_list)){
  calls$mut_id[i] <- call_list[i]
  calls$sample[i] <- strsplit(calls$mut_id[i], split = "_")[[1]][[1]]
  calls$pos[i] <- as.numeric(strsplit(calls$mut_id[i], split = "_")[[1]][[2]])
  calls$ref[i] <- strsplit(calls$mut_id[i], split = "_")[[1]][[3]]
  calls$mut[i] <- strsplit(calls$mut_id[i], split = "_")[[1]][[4]]
  if(nchar(calls$ref[i]) == 1){
    calls$r1_xfw[i] <- run_1_coverage[which(sample_list == calls$sample[i]),
                                            calls$pos[i],
                                            which(ref_list == calls$mut[i])]
    calls$r1_xbw[i] <- run_1_coverage[which(sample_list == calls$sample[i]),
                                            calls$pos[i],
                                            (which(ref_list == calls$mut[i]) + 6)]
    calls$r1_nfw[i] <- sum(run_1_coverage[which(sample_list == calls$sample[i]),
                                                calls$pos[i],
                                                1:6])
    calls$r1_nbw[i] <- sum(run_1_coverage[which(sample_list == calls$sample[i]),
                                                calls$pos[i],
                                                7:12])
    calls$r1_vaf[i] <- signif((as.numeric(calls$r1_xfw[i]) + as.numeric(calls$r1_xbw[i])) / (as.numeric(calls$r1_nfw[i]) + as.numeric(calls$r1_nbw[i])), digits = 2)
    
    calls$r2_xfw[i] <- run_2_coverage[which(sample_list == calls$sample[i]),
                                            calls$pos[i],
                                            which(ref_list == calls$mut[i])]
    calls$r2_xbw[i] <- run_2_coverage[which(sample_list == calls$sample[i]),
                                            calls$pos[i],
                                            (which(ref_list == calls$mut[i]) + 6)]
    calls$r2_nfw[i] <- sum(run_2_coverage[which(sample_list == calls$sample[i]),
                                                calls$pos[i],
                                                1:6])
    calls$r2_nbw[i] <- sum(run_2_coverage[which(sample_list == calls$sample[i]),
                                                calls$pos[i],
                                                7:12])
    calls$r2_vaf[i] <- signif((as.numeric(calls$r2_xfw[i]) + as.numeric(calls$r2_xbw[i])) / (as.numeric(calls$r2_nfw[i]) + as.numeric(calls$r2_nbw[i])),digits = 2)
    
    row_num_1 <- which(initial_sw_1_1$mut_id == calls$mut_id[i])
    if(length(row_num_1) > 0){
      calls$r1_pval[i] <- signif(initial_sw_1_1$pval[row_num_1], digits = 2)
    }
    row_num_2 <- which(initial_sw_1_2$mut_id == calls$mut_id[i])
    if(length(row_num_2) > 0){
      calls$r2_pval[i] <- signif(initial_sw_1_2$pval[row_num_2], digits = 2)
    }
  } else{
    change_length <- nchar(calls$ref[i])
    calls$r1_xfw[i] <- paste(run_1_coverage[which(sample_list == calls$sample[i]),
                                                  calls$pos[i]:(calls$pos[i]+change_length-1),
                                                  which(ref_list == calls$mut[i])],collapse = ";")
    calls$r1_xbw[i] <- paste(run_1_coverage[which(sample_list == calls$sample[i]),
                                                  calls$pos[i]:(calls$pos[i]+change_length-1),
                                                  (which(ref_list == calls$mut[i]) + 6)],collapse = ";")
    calls$r1_nfw[i] <- paste(rowSums(run_1_coverage[which(sample_list == calls$sample[i]),
                                                          calls$pos[i]:(calls$pos[i]+change_length-1),
                                                          1:6]),collapse = ";")
    calls$r1_nbw[i] <- paste(rowSums(run_1_coverage[which(sample_list == calls$sample[i]),
                                                          calls$pos[i]:(calls$pos[i]+change_length-1),
                                                          7:12]),collapse = ";")
    calls$r1_vaf[i] <- paste(signif((run_1_coverage[which(sample_list == calls$sample[i]),calls$pos[i]:(calls$pos[i]+change_length-1),which(ref_list == calls$mut[i])]+run_1_coverage[which(sample_list == calls$sample[i]),calls$pos[i]:(calls$pos[i]+change_length-1),(which(ref_list == calls$mut[i]) + 6)])/(rowSums(run_1_coverage[which(sample_list == calls$sample[i]),calls$pos[i]:(calls$pos[i]+change_length-1),1:6])+rowSums(run_1_coverage[which(sample_list == calls$sample[i]),calls$pos[i]:(calls$pos[i]+change_length-1),7:12])), digits = 2),collapse = ";")
    
    calls$r2_xfw[i] <- paste(run_2_coverage[which(sample_list == calls$sample[i]),
                                                  calls$pos[i]:(calls$pos[i]+change_length-1),
                                                  which(ref_list == calls$mut[i])],collapse = ";")
    calls$r2_xbw[i] <- paste(run_2_coverage[which(sample_list == calls$sample[i]),
                                                  calls$pos[i]:(calls$pos[i]+change_length-1),
                                                  (which(ref_list == calls$mut[i]) + 6)],collapse = ";")
    calls$r2_nfw[i] <- paste(rowSums(run_2_coverage[which(sample_list == calls$sample[i]),
                                                          calls$pos[i]:(calls$pos[i]+change_length-1),
                                                          1:6]),collapse = ";")
    calls$r2_nbw[i] <- paste(rowSums(run_2_coverage[which(sample_list == calls$sample[i]),
                                                          calls$pos[i]:(calls$pos[i]+change_length-1),
                                                          7:12]),collapse = ";")
    calls$r2_vaf[i] <- paste(signif((run_2_coverage[which(sample_list == calls$sample[i]),calls$pos[i]:(calls$pos[i]+change_length-1),which(ref_list == calls$mut[i])]+run_2_coverage[which(sample_list == calls$sample[i]),calls$pos[i]:(calls$pos[i]+change_length-1),(which(ref_list == calls$mut[i]) + 6)])/(rowSums(run_2_coverage[which(sample_list == calls$sample[i]),calls$pos[i]:(calls$pos[i]+change_length-1),1:6])+rowSums(run_2_coverage[which(sample_list == calls$sample[i]),calls$pos[i]:(calls$pos[i]+change_length-1),7:12])), digits = 2),collapse = ";")
    
    pval_1_array <- vector(length = change_length)
    pval_1_array[] <- NA
    
    pval_2_array <- vector(length = change_length)
    pval_2_array[] <- NA
    
    for(j in 1:change_length){
      if(nchar(calls$mut[i]) == nchar(calls$ref[i])){
        temp_row_num_1 <- which(initial_sw_1_1$sampleID == calls$sample[i] &
                                  initial_sw_1_1$pos == (calls$pos[i] + j - 1) &
                                  initial_sw_1_1$mut == substr(x = calls$mut[i],start = j, stop = j))
      }else{
        temp_row_num_1 <- which(initial_sw_1_1$sampleID == calls$sample[i] &
                                  initial_sw_1_1$pos == (calls$pos[i] + j - 1) &
                                  initial_sw_1_1$mut == calls$mut[i])
        
      }
           if(length(temp_row_num_1) > 0){
        pval_1_array[j] <- signif(initial_sw_1_1$pval[temp_row_num_1],digits = 2)
      }
      
      if(nchar(calls$mut[i]) == nchar(calls$ref[i])){
        temp_row_num_2 <- which(initial_sw_1_2$sampleID == calls$sample[i] &
                                  initial_sw_1_2$pos == (calls$pos[i] + j - 1) &
                                  initial_sw_1_2$mut == substr(x = calls$mut[i],start = j, stop = j))
      }else{
        temp_row_num_2 <- which(initial_sw_1_2$sampleID == calls$sample[i] &
                                  initial_sw_1_2$pos == (calls$pos[i] + j - 1) &
                                  initial_sw_1_2$mut == calls$mut[i])
        
      }
      if(length(temp_row_num_2) > 0){
        pval_2_array[j] <- signif(initial_sw_1_2$pval[temp_row_num_2], digits = 2)
      }
    }
    calls$r1_pval[i] <- paste(pval_1_array, collapse = ";")
    calls$r2_pval[i] <- paste(pval_2_array, collapse = ";")
  }
  
  row_num_strict_1 <- which(sw_1_1$mut_id == calls$mut_id[i])
  if(length(row_num_strict_1) > 0){
    calls$r1_strict[i] <- T
    calls$r1_qval[i] <- signif(sw_1_1$qval[row_num_strict_1],2)
  } else{
    calls$r1_strict[i] <- F
  }
  row_num_strict_2 <- which(sw_1_2$mut_id == calls$mut_id[i])
  if(length(row_num_strict_2) > 0){
    calls$r2_strict[i] <- T
    calls$r2_qval[i] <-  signif(sw_1_2$qval[row_num_strict_2],2)
  } else{
    calls$r2_strict[i] <- F
  }
}

multinuc_subs <- calls[which(nchar(calls$ref) == nchar(calls$mut) & nchar(calls$ref) != 1),]
other_calls <- calls[which(!(nchar(calls$ref) == nchar(calls$mut) & nchar(calls$ref) != 1)),]

for(i in 1:nrow(multinuc_subs)){
  change_length <- nchar(multinuc_subs$ref[i])
  
  temp_r1_xfw <- vector()
  temp_r1_xbw <- vector()
  temp_r1_nfw <- vector()
  temp_r1_nbw <- vector()
  temp_r1_vaf <- vector()
  
  temp_r2_xfw <- vector()
  temp_r2_xbw <- vector()
  temp_r2_nfw <- vector()
  temp_r2_nbw <- vector()
  temp_r2_vaf <- vector()
  
  for(j in 1:change_length){
    temp_r1_xfw <- c(temp_r1_xfw, run_1_coverage[which(sample_list == multinuc_subs$sample[i]),
                                                 (multinuc_subs$pos[i]+j-1),
                                                 which(ref_list == substr(multinuc_subs$mut[i],j,j))])  
  
    temp_r1_xbw <- c(temp_r1_xbw, run_1_coverage[which(sample_list == multinuc_subs$sample[i]),
                                                 (multinuc_subs$pos[i]+j-1),
                                                 (which(ref_list == substr(multinuc_subs$mut[i],j,j)) + 6)])
    
    temp_r1_nfw <- c(temp_r1_nfw, sum(run_1_coverage[which(sample_list == multinuc_subs$sample[i]),
                                                 (multinuc_subs$pos[i]+j-1),
                                                 1:6]))
    temp_r1_nbw <- c(temp_r1_nbw, sum(run_1_coverage[which(sample_list == multinuc_subs$sample[i]),
                                                     (multinuc_subs$pos[i]+j-1),
                                                     7:12]))
    temp_r1_vaf <- c(temp_r1_vaf, signif((temp_r1_xfw[j] + temp_r1_xbw[j]) / (temp_r1_nfw[j] + temp_r1_nbw[j]), digits = 2))
      
    temp_r2_xfw <- c(temp_r2_xfw, run_2_coverage[which(sample_list == multinuc_subs$sample[i]),
                                                 (multinuc_subs$pos[i]+j-1),
                                                 which(ref_list == substr(multinuc_subs$mut[i],j,j))])  
    
    temp_r2_xbw <- c(temp_r2_xbw, run_2_coverage[which(sample_list == multinuc_subs$sample[i]),
                                                 (multinuc_subs$pos[i]+j-1),
                                                 (which(ref_list == substr(multinuc_subs$mut[i],j,j)) + 6)])
    temp_r2_nfw <- c(temp_r2_nfw, sum(run_2_coverage[which(sample_list == multinuc_subs$sample[i]),
                                                     (multinuc_subs$pos[i]+j-1),
                                                     1:6]))
    temp_r2_nbw <- c(temp_r2_nbw, sum(run_2_coverage[which(sample_list == multinuc_subs$sample[i]),
                                                     (multinuc_subs$pos[i]+j-1),
                                                     7:12]))
    temp_r2_vaf <- c(temp_r2_vaf, signif((temp_r2_xfw[j] + temp_r2_xbw[j]) / (temp_r2_nfw[j] + temp_r2_nbw[j]), digits = 2))
  }
  multinuc_subs$r1_xfw[i] <- paste(temp_r1_xfw, collapse = ";")
  multinuc_subs$r1_xbw[i] <- paste(temp_r1_xbw, collapse = ";")
  multinuc_subs$r1_vaf[i] <- paste(temp_r1_vaf, collapse = ";")
  multinuc_subs$r2_xfw[i] <- paste(temp_r2_xfw, collapse = ";")
  multinuc_subs$r2_xbw[i] <- paste(temp_r2_xbw, collapse = ";")
  multinuc_subs$r2_vaf[i] <- paste(temp_r2_vaf, collapse = ";")
}

calls <- rbind(other_calls,multinuc_subs)
calls <- calls[order(calls$sample,calls$pos),]

saveRDS(calls, file = "shearwater_intermediate_intersect_muts_both_runs_annotated_deltaVAF0.05merge.rds")

singlenuc_calls <- calls[which(nchar(calls$ref) == 1),]
singlenuc_calls$r1_xfw <- as.numeric(singlenuc_calls$r1_xfw)
singlenuc_calls$r1_xbw <- as.numeric(singlenuc_calls$r1_xbw)
singlenuc_calls$r1_nfw <- as.numeric(singlenuc_calls$r1_nfw)
singlenuc_calls$r1_nbw <- as.numeric(singlenuc_calls$r1_nbw)
singlenuc_calls$r1_vaf <- as.numeric(singlenuc_calls$r1_vaf)

singlenuc_calls$r2_xfw <- as.numeric(singlenuc_calls$r2_xfw)
singlenuc_calls$r2_xbw <- as.numeric(singlenuc_calls$r2_xbw)
singlenuc_calls$r2_nfw <- as.numeric(singlenuc_calls$r2_nfw)
singlenuc_calls$r2_nbw <- as.numeric(singlenuc_calls$r2_nbw)
singlenuc_calls$r2_vaf <- as.numeric(singlenuc_calls$r2_vaf)

singlenuc_calls$both_vaf <- (singlenuc_calls$r1_xfw + singlenuc_calls$r1_xbw + singlenuc_calls$r2_xfw + singlenuc_calls$r2_xbw) / (singlenuc_calls$r1_nfw + singlenuc_calls$r1_nbw + singlenuc_calls$r2_nfw + singlenuc_calls$r2_nbw)

multinuc_calls <- calls[which(nchar(calls$ref) > 1),]
multinuc_calls_median <- multinuc_calls

for(i in 1:nrow(multinuc_calls_median)){
  multinuc_calls_median$r1_xfw[i] <- floor(median(as.numeric(unlist(strsplit(multinuc_calls$r1_xfw[i], split = ";")))))
  multinuc_calls_median$r1_xbw[i] <- floor(median(as.numeric(unlist(strsplit(multinuc_calls$r1_xbw[i], split = ";")))))
  multinuc_calls_median$r1_nfw[i] <- ceiling(median(as.numeric(unlist(strsplit(multinuc_calls$r1_nfw[i], split = ";")))))
  multinuc_calls_median$r1_nbw[i] <- ceiling(median(as.numeric(unlist(strsplit(multinuc_calls$r1_nbw[i], split = ";")))))
  multinuc_calls_median$r1_vaf[i] <- mean(as.numeric(unlist(strsplit(multinuc_calls$r1_vaf[i], split = ";"))))
  
  multinuc_calls_median$r2_xfw[i] <- floor(median(as.numeric(unlist(strsplit(multinuc_calls$r2_xfw[i], split = ";")))))
  multinuc_calls_median$r2_xbw[i] <- floor(median(as.numeric(unlist(strsplit(multinuc_calls$r2_xbw[i], split = ";")))))
  multinuc_calls_median$r2_nfw[i] <- ceiling(median(as.numeric(unlist(strsplit(multinuc_calls$r2_nfw[i], split = ";")))))
  multinuc_calls_median$r2_nbw[i] <- ceiling(median(as.numeric(unlist(strsplit(multinuc_calls$r2_nbw[i], split = ";")))))
  multinuc_calls_median$r2_vaf[i] <- mean(as.numeric(unlist(strsplit(multinuc_calls$r2_vaf[i], split = ";"))))
}

multinuc_calls_median$r1_xfw <- as.numeric(multinuc_calls_median$r1_xfw)
multinuc_calls_median$r1_xbw <- as.numeric(multinuc_calls_median$r1_xbw)
multinuc_calls_median$r1_nfw <- as.numeric(multinuc_calls_median$r1_nfw)
multinuc_calls_median$r1_nbw <- as.numeric(multinuc_calls_median$r1_nbw)
multinuc_calls_median$r1_vaf <- as.numeric(multinuc_calls_median$r1_vaf)

multinuc_calls_median$r2_xfw <- as.numeric(multinuc_calls_median$r2_xfw)
multinuc_calls_median$r2_xbw <- as.numeric(multinuc_calls_median$r2_xbw)
multinuc_calls_median$r2_nfw <- as.numeric(multinuc_calls_median$r2_nfw)
multinuc_calls_median$r2_nbw <- as.numeric(multinuc_calls_median$r2_nbw)
multinuc_calls_median$r2_vaf <- as.numeric(multinuc_calls_median$r2_vaf)

multinuc_calls_median$both_vaf <- (multinuc_calls_median$r1_xfw + multinuc_calls_median$r1_xbw + multinuc_calls_median$r2_xfw + multinuc_calls_median$r2_xbw) / (multinuc_calls_median$r1_nfw + multinuc_calls_median$r1_nbw + multinuc_calls_median$r2_nfw + multinuc_calls_median$r2_nbw)

calls_median <- rbind(singlenuc_calls,multinuc_calls_median)
calls_median <- calls_median[order(calls_median$sample,calls_median$pos),]

saveRDS(calls_median, file = "shearwater_intermediate_intersect_muts_both_runs_annotated_median_values_deltaVAF0.05_merge.rds")
