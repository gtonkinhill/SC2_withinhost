library("Rsamtools")

sw_1_1 <- readRDS("shearwater_muts_all_samples_qval0.05_mincov100x_deltaVAF0.05merge_panel_1_run_1.rds")
sw_1_2 <- readRDS("shearwater_muts_all_samples_qval0.05_mincov100x_deltaVAF0.05merge_panel_1_run_2.rds")

initial_sw_1_1 <- readRDS("initial_shearwater_calls_background_panel_1_run_1.rds")
initial_sw_1_2 <- readRDS("initial_shearwater_calls_background_panel_1_run_2.rds")

pval_0.01_sw_1_1 <- initial_sw_1_1[which(initial_sw_1_1$pval <= 0.01),]
pval_0.01_sw_1_2 <- initial_sw_1_2[which(initial_sw_1_2$pval <= 0.01),]

pval_0.01_sw_1_1$mut_site <- paste(pval_0.01_sw_1_1$pos,pval_0.01_sw_1_1$ref,pval_0.01_sw_1_1$mut,sep = "_")
pval_0.01_sw_1_1$mut_id <- paste(pval_0.01_sw_1_1$sampleID,pval_0.01_sw_1_1$mut_site,sep = "_")
pval_0.01_sw_1_2$mut_site <- paste(pval_0.01_sw_1_2$pos,pval_0.01_sw_1_2$ref,pval_0.01_sw_1_2$mut,sep = "_")
pval_0.01_sw_1_2$mut_id <- paste(pval_0.01_sw_1_2$sampleID,pval_0.01_sw_1_2$mut_site,sep = "_")

mutreads_5_sw_1_1 <- initial_sw_1_1[which((initial_sw_1_1$xfw + initial_sw_1_1$xbw) >= 5),]
mutreads_5_sw_1_2 <- initial_sw_1_2[which((initial_sw_1_2$xfw + initial_sw_1_2$xbw) >= 5),]

mutreads_5_sw_1_1$mut_site <- paste(mutreads_5_sw_1_1$pos,mutreads_5_sw_1_1$ref,mutreads_5_sw_1_1$mut,sep = "_")
mutreads_5_sw_1_1$mut_id <- paste(mutreads_5_sw_1_1$sampleID,mutreads_5_sw_1_1$mut_site,sep = "_")
mutreads_5_sw_1_2$mut_site <- paste(mutreads_5_sw_1_2$pos,mutreads_5_sw_1_2$ref,mutreads_5_sw_1_2$mut,sep = "_")
mutreads_5_sw_1_2$mut_id <- paste(mutreads_5_sw_1_2$sampleID,mutreads_5_sw_1_2$mut_site,sep = "_")

genomeFile <- ("MN908947.3.fasta")

trinuc_count_minus <- readRDS("MN908947.3_trinuc_count_py_on_minus_strand.rds")
trinuc_count_plus <- readRDS("MN908947.3_trinuc_count_py_on_plus_strand.rds")

trinucleotide_plot = function (mutations, file_name, analysis_type, analysis_region) {
  mutations <- unique(mutations[,c("chr","pos","ref","mut")])
  mutations$ref <- as.character(mutations$ref)
  mutations$mut <- as.character(mutations$mut)
  mutations <- mutations[(mutations$ref %in% c("A","C","G","T")) & (mutations$mut %in% c("A","C","G","T")),]
  
  if(analysis_region == "exc_regions"){
    mutations <- mutations[which(mutations$pos %in% c(1:54,4996:5259,6847:7058,7093:7332,11669:11889,16486:16770,19276:19570,19912:20496,21147:21386,22325:22542,27512:27808,29837:29903)),]
  }else if(analysis_region == "inc_regions"){
    mutations <- mutations[which(!(mutations$pos %in% c(1:54,4996:5259,6847:7058,7093:7332,11669:11889,16486:16770,19276:19570,19912:20496,21147:21386,22325:22542,27512:27808,29837:29903))),]
  }else if(analysis_region != "whole_genome"){
    mutations <- NULL
  }
  
  if(nrow(mutations) > 0){
    mutations$trinuc_ref = as.vector(scanFa(genomeFile, GRanges(mutations$chr, IRanges(mutations$pos-1, mutations$pos+1))))
    
    # 2. Annotating the mutation from the pyrimidine base
    ntcomp = c(T="A",G="C",C="G",A="T")
    mutations$sub = paste(mutations$ref,mutations$mut,sep=">")
    mutations$trinuc_ref_py = mutations$trinuc_ref
    for (j in 1:nrow(mutations)) {
      if (mutations$ref[j] %in% c("A","G")) { # Purine base
        mutations$sub[j] = paste(ntcomp[mutations$ref[j]],ntcomp[mutations$mut[j]],sep=">")
        mutations$trinuc_ref_py[j] = paste(ntcomp[rev(strsplit(mutations$trinuc_ref[j],split="")[[1]])],collapse="")
      }
    }
    
    # 3. Counting subs
    freqs_minus = table(paste(mutations$sub[which(mutations$ref %in% c("A","G"))],paste(substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("A","G"))],1,1),substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("A","G"))],3,3),sep="-"),sep=","))
    freqs_plus = table(paste(mutations$sub[which(mutations$ref %in% c("C","T"))],paste(substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("C","T"))],1,1),substr(mutations$trinuc_ref_py[which(mutations$ref %in% c("C","T"))],3,3),sep="-"),sep=","))
    
    sub_vec = c("C>A","C>G","C>T","T>A","T>C","T>G")
    ctx_vec = paste(rep(c("A","C","G","T"),each=4),rep(c("A","C","G","T"),times=4),sep="-")
    full_vec = paste(rep(sub_vec,each=16),rep(ctx_vec,times=6),sep=",")
    freqs_minus_full = freqs_minus[full_vec]; freqs_minus_full[is.na(freqs_minus_full)] = 0; names(freqs_minus_full) = full_vec
    freqs_plus_full = freqs_plus[full_vec]; freqs_plus_full[is.na(freqs_plus_full)] = 0; names(freqs_plus_full) = full_vec
    
    if(analysis_type == "obs_exp"){
      if(analysis_region == "exc_regions"){
        minus_base_freqs <- (trinuc_count_minus[which(trinuc_count_minus$name == "exc_regions_only"),6:37] / sum(trinuc_count_minus[which(trinuc_count_minus$name == "exc_regions_only"),6:37]))
        exp_minus_counts <- sum(freqs_minus) * as.numeric(c(rep(minus_base_freqs[1:16]/3,times = 3),rep(minus_base_freqs[17:32]/3,times = 3)))
        freqs_minus_full <- freqs_minus_full / exp_minus_counts
        
        plus_base_freqs <- (trinuc_count_plus[which(trinuc_count_plus$name == "exc_regions_only"),6:37] / sum(trinuc_count_plus[which(trinuc_count_plus$name == "exc_regions_only"),6:37]))
        exp_plus_counts <- sum(freqs_plus) * as.numeric(c(rep(plus_base_freqs[1:16]/3,times = 3),rep(plus_base_freqs[17:32]/3,times = 3)))
        freqs_plus_full <- freqs_plus_full / exp_plus_counts
      }else if(analysis_region == "inc_regions"){
        minus_base_freqs <- (trinuc_count_minus[which(trinuc_count_minus$name == "inc_regions_only"),6:37] / sum(trinuc_count_minus[which(trinuc_count_minus$name == "inc_regions_only"),6:37]))
        exp_minus_counts <- sum(freqs_minus) * as.numeric(c(rep(minus_base_freqs[1:16]/3,times = 3),rep(minus_base_freqs[17:32]/3,times = 3)))
        freqs_minus_full <- freqs_minus_full / exp_minus_counts
        
        plus_base_freqs <- (trinuc_count_plus[which(trinuc_count_plus$name == "inc_regions_only"),6:37] / sum(trinuc_count_plus[which(trinuc_count_plus$name == "inc_regions_only"),6:37]))
        exp_plus_counts <- sum(freqs_plus) * as.numeric(c(rep(plus_base_freqs[1:16]/3,times = 3),rep(plus_base_freqs[17:32]/3,times = 3)))
        freqs_plus_full <- freqs_plus_full / exp_plus_counts
      }else if(analysis_region == "whole_genome"){
        minus_base_freqs <- (trinuc_count_minus[which(trinuc_count_minus$name == "whole_genome"),6:37] / sum(trinuc_count_minus[which(trinuc_count_minus$name == "whole_genome"),6:37]))
        exp_minus_counts <- sum(freqs_minus) * as.numeric(c(rep(minus_base_freqs[1:16]/3,times = 3),rep(minus_base_freqs[17:32]/3,times = 3)))
        freqs_minus_full <- freqs_minus_full / exp_minus_counts
        
        plus_base_freqs <- (trinuc_count_plus[which(trinuc_count_plus$name == "whole_genome"),6:37] / sum(trinuc_count_plus[which(trinuc_count_plus$name == "whole_genome"),6:37]))
        exp_plus_counts <- sum(freqs_plus) * as.numeric(c(rep(plus_base_freqs[1:16]/3,times = 3),rep(plus_base_freqs[17:32]/3,times = 3)))
        freqs_plus_full <- freqs_plus_full / exp_plus_counts
      }else{
        freqs_minus_full = NULL
        freqs_plus_full = NULL
      }
    }
    
    xstr = paste(substr(full_vec,5,5), substr(full_vec,1,1), substr(full_vec,7,7), sep="")
    
    dev.new(width=10,height=4)
    colvec = rep(c("dodgerblue","black","red","grey70","olivedrab3","plum2"),each=16)
    y_minus = freqs_minus_full; y_plus = freqs_plus_full;
    
    maxy = max(c(y_minus,y_plus), na.rm = T)
    
    if(analysis_type == "obs_exp"){
      ylab = "Mutation frequency (Obs/Exp)"
    }else{
      ylab = "Mutation count"
    }
    
    h_plus = barplot(y_plus, las=2, col=colvec, border=NA, ylim=c(-maxy*1.5,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab=ylab)
    h_minus = barplot(-y_minus, las=2, col=colvec, border=NA, ylim=c(-maxy*1.5,maxy*1.5), space=1, cex.names=0.6, names.arg=xstr, ylab=ylab, add = T)
    
    segments(y0 = maxy*1.5, y1 = maxy*1.5, x0 = 0.5, x1 = 192.5,  col = "black")
    segments(y0 = -maxy*1.5, y1 = -maxy*1.5, x0 = 0.5, x1 = 192.5,  col = "black")
    segments(y0 = 0, y1 = 0, x0 = 0.5, x1 = 192.5,  col = "black")
    abline(v = 0.5, col = "black")
    abline(v = 32.5, col = "black")
    abline(v = 64.5, col = "black")
    abline(v = 96.5, col = "black")
    abline(v = 128.5, col = "black")
    abline(v = 160.5, col = "black")
    abline(v = 192.5, col = "black")
    
    
    for (j in 1:length(sub_vec)) {
      xpos = h_minus[c((j-1)*16+1,j*16)]
      rect(xpos[1]-0.5, maxy*1.25, xpos[2]+0.5, maxy*1.15, border=NA, col=colvec[j*16])
      text(x=mean(xpos), y=maxy*1.15, pos=3, labels=sub_vec[j])
    }    
    dev.copy(pdf,file_name,width=12,height=5)
    dev.off()
    dev.off()
  }
}

lower_bound <- c(0,0,0.01,0.02,0.05,0.1,0.25,0.5)
upper_bound <- c(1,0.01,0.02,0.05,0.1,0.25,0.5,1)

union_muts <- rbind(sw_1_1,sw_1_2)
union_muts <- union_muts[order(union_muts$vaf, decreasing = T),]
union_muts <- union_muts[which(!(duplicated(union_muts$mut_id))),]

union_mut_list <- union_muts$mut_id
strict_intersect_mut_list <- union_mut_list[which(union_mut_list %in% sw_1_1$mut_id &
                                             union_mut_list %in% sw_1_2$mut_id)]

strict_intersect_muts <- rbind(sw_1_1,sw_1_2)
strict_intersect_muts <- strict_intersect_muts[which(strict_intersect_muts$mut_id %in% strict_intersect_mut_list),]
strict_intersect_muts <- strict_intersect_muts[order(strict_intersect_muts$vaf, decreasing = T),]
strict_intersect_muts <- strict_intersect_muts[which(!(duplicated(strict_intersect_muts$mut_id))),]


intermediate_intersect_mut_list <- union_mut_list[which((union_mut_list %in% sw_1_1$mut_id &
                                                      union_mut_list %in% pval_0.01_sw_1_2$mut_id) |
                                                     (union_mut_list %in% sw_1_2$mut_id &
                                                        union_mut_list %in% pval_0.01_sw_1_1$mut_id))]

multinuc_union_mut_list <- union_mut_list[which(nchar(unlist(strsplit(union_mut_list, split = "_"))[seq(from = 3, by = 4, to = (4*length(union_mut_list) - 1))]) > 1)]

multinuc_intermediate_muts <- vector()
multinuc_relaxed_muts <- vector()

for(i in 23348:length(multinuc_union_mut_list)){
  temp_sample <- unlist(strsplit(multinuc_union_mut_list[i], split = "_"))[1]
  temp_pos <- as.numeric(unlist(strsplit(multinuc_union_mut_list[i], split = "_"))[2])
  temp_ref <- unlist(strsplit(multinuc_union_mut_list[i], split = "_"))[3]
  temp_mut <- unlist(strsplit(multinuc_union_mut_list[i], split = "_"))[4]
  change_length <- nchar(temp_ref)
  
  qval_array <- vector(length = 2)
  if(multinuc_union_mut_list[i] %in% sw_1_1$mut_id){
    qval_array[1] <- T
  }
  if(multinuc_union_mut_list[i] %in% sw_1_2$mut_id){
    qval_array[2] <- T
  }
  
  if(all(qval_array)){
    multinuc_intermediate_muts <- c(multinuc_intermediate_muts,multinuc_union_mut_list[i])
    multinuc_relaxed_muts <- c(multinuc_relaxed_muts, multinuc_union_mut_list[i])
  } else{
    pval_array <- array(data = F, dim = c(2,change_length))
    dep_array <- array(data = F, dim = c(2, change_length))
    for(j in 1:change_length){
      if(nchar(temp_ref) != nchar(temp_mut)){
        pval_array[1,j] <- paste(temp_sample,(temp_pos + j - 1),substr(temp_ref,j,j),temp_mut,sep = "_")%in% pval_0.01_sw_1_1$mut_id
        pval_array[2,j] <- paste(temp_sample,(temp_pos + j - 1),substr(temp_ref,j,j),temp_mut,sep = "_")%in% pval_0.01_sw_1_2$mut_id
        dep_array[1,j] <- paste(temp_sample,(temp_pos + j - 1),substr(temp_ref,j,j),temp_mut,sep = "_")%in% mutreads_5_sw_1_1$mut_id
        dep_array[2,j] <- paste(temp_sample,(temp_pos + j - 1),substr(temp_ref,j,j),temp_mut,sep = "_")%in% mutreads_5_sw_1_2$mut_id
      } else{
        pval_array[1,j] <- paste(temp_sample,(temp_pos + j - 1),substr(temp_ref,j,j),substr(temp_mut,j,j),sep = "_")%in% pval_0.01_sw_1_1$mut_id
        pval_array[2,j] <- paste(temp_sample,(temp_pos + j - 1),substr(temp_ref,j,j),substr(temp_mut,j,j),sep = "_")%in% pval_0.01_sw_1_2$mut_id
        dep_array[1,j] <- paste(temp_sample,(temp_pos + j - 1),substr(temp_ref,j,j),substr(temp_mut,j,j),sep = "_")%in% mutreads_5_sw_1_1$mut_id
        dep_array[2,j] <- paste(temp_sample,(temp_pos + j - 1),substr(temp_ref,j,j),substr(temp_mut,j,j),sep = "_")%in% mutreads_5_sw_1_2$mut_id
      }
    }
    if((qval_array[1] & any(pval_array[2,])) | (qval_array[2] & any(pval_array[1,]))){
      multinuc_intermediate_muts <- c(multinuc_intermediate_muts,multinuc_union_mut_list[i])
    }
    if((qval_array[1] & any(dep_array[2,])) | (qval_array[2] & any(dep_array[1,]))){
      multinuc_relaxed_muts <- c(multinuc_relaxed_muts,multinuc_union_mut_list[i])
    }
  }
}

intermediate_intersect_muts <- rbind(sw_1_1,sw_1_2)
intermediate_intersect_muts <- intermediate_intersect_muts[which(intermediate_intersect_muts$mut_id %in% intermediate_intersect_mut_list | intermediate_intersect_muts$mut_id %in% multinuc_intermediate_muts),]
intermediate_intersect_muts <- intermediate_intersect_muts[order(intermediate_intersect_muts$vaf, decreasing = T),]
intermediate_intersect_muts <- intermediate_intersect_muts[which(!(duplicated(intermediate_intersect_muts$mut_id))),]
saveRDS(intermediate_intersect_muts,"shearwater_intermediate_intersect_muts_deltaVAF0.05merge.rds")
saveRDS(c(intermediate_intersect_mut_list,multinuc_intermediate_muts),"shearwater_intermediate_intersect_mut_list_deltaVAF0.05merge.rds")

relaxed_intersect_mut_list <- union_mut_list[which((union_mut_list %in% sw_1_1$mut_id &
                                                      union_mut_list %in% mutreads_5_sw_1_2$mut_id) |
                                                     (union_mut_list %in% sw_1_2$mut_id &
                                                        union_mut_list %in% mutreads_5_sw_1_1$mut_id))]

relaxed_intersect_muts <- rbind(sw_1_1,sw_1_2)
relaxed_intersect_muts <- relaxed_intersect_muts[which(relaxed_intersect_muts$mut_id %in% relaxed_intersect_mut_list | relaxed_intersect_muts$mut_id %in% multinuc_relaxed_muts),]
relaxed_intersect_muts <- relaxed_intersect_muts[order(relaxed_intersect_muts$vaf, decreasing = T),]
relaxed_intersect_muts <- relaxed_intersect_muts[which(!(duplicated(relaxed_intersect_muts$mut_id))),]
saveRDS(relaxed_intersect_muts,"shearwater_relaxed_intersect_muts_deltaVAF0.05merge.rds")
saveRDS(c(relaxed_intersect_mut_list,multinuc_relaxed_muts),"shearwater_relaxed_intersect_mut_list_deltaVAF0.05merge.rds")

dataset_list <- c("strict_intersect_muts","intermediate_intersect_muts","relaxed_intersect_muts","union_muts")

for(j in 1:length(dataset_list)){
  dataset <- eval(parse(text = dataset_list[j]))
  
  for(i in 1:length(lower_bound)){
    if(upper_bound[i] == 1){
      mutations <- dataset[which(dataset$vaf >= lower_bound[i] & dataset$vaf <= upper_bound[i]),]
    } else{
      mutations <- dataset[which(dataset$vaf >= lower_bound[i] & dataset$vaf < upper_bound[i]),]
    }
    trinucleotide_plot(mutations,paste("shearwater_all_samples_qval0.05_mindep100x_",dataset_list[j],"_whole_genome_highest_vaf_from",lower_bound[i],"to",upper_bound[i],"_trinucleotide_counts.pdf",sep=""),"counts","whole_genome")
    #    trinucleotide_plot(mutations,paste("shearwater_",dataset_list[j],"_inc_regions_highest_vaf_from",lower_bound[i],"to",upper_bound[i],"_trinucleotide_counts.pdf",sep=""),"counts","inc_regions")
    #    trinucleotide_plot(mutations,paste("shearwater_",dataset_list[j],"_exc_regions_highest_vaf_from",lower_bound[i],"to",upper_bound[i],"_trinucleotide_counts.pdf",sep=""),"counts","exc_regions")
    
    trinucleotide_plot(mutations,paste("shearwater_all_samples_qval0.05_mindep100x_",dataset_list[j],"_whole_genome_highest_vaf_from",lower_bound[i],"to",upper_bound[i],"_trinucleotide_obs_exp.pdf",sep=""),"obs_exp","whole_genome")
    #    trinucleotide_plot(mutations,paste("shearwater_",dataset_list[j],"_inc_regions_highest_vaf_from",lower_bound[i],"to",upper_bound[i],"_trinucleotide_obs_exp.pdf",sep=""),"obs_exp","inc_regions")
    #    trinucleotide_plot(mutations,paste("shearwater_",dataset_list[j],"_exc_regions_highest_vaf_from",lower_bound[i],"to",upper_bound[i],"_trinucleotide_obs_exp.pdf",sep=""),"obs_exp","exc_regions")
  }
}  
