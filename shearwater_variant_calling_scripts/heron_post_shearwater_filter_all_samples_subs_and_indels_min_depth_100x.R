shearwater <- readRDS("initial_shearwater_calls_background_panel_1_run_1.rds")

deltavaf <- 0.05 # Maximum difference in VAF allowed to merge mutations

sample_list <- read.table("sample_list_1181.tsv", sep = "\t", header = F, stringsAsFactors = F)
sample_list <- sample_list$V1

coverage <- readRDS("run_1_coverage.rds")

mutations <- shearwater[which((shearwater$nfw + shearwater$nbw) >= 100 & shearwater$pval <= 0.001),]
mutations$sampleID <- as.character(mutations$sampleID)

sites_with_100x <- 0
for(i in 1:length(sample_list)){
  sites_with_100x <- sites_with_100x + length(which(rowSums(coverage[i,,]) >= 100))
}

mutations$vaf = (mutations$xfw+mutations$xbw)/(mutations$nfw+mutations$nbw)

indels         <- mutations[which(mutations$mut=="-"),];
indels         <- indels[order(indels$sampleID, indels$chr, indels$pos),];
i              <- 1;
group_counter  <- 1;
while(i <= nrow(indels)) {
  indel_from  <- indels[i,"pos"];
  indel_to    <- indel_from;
  from_index  <- i;
  to_index    <- i;
  deleted_seq <- indels[i,"ref"];
  sample_f    <- indels[i,"sampleID"];
  i <- i+1;
  for(j in c(i:nrow(indels))) {
    if(j>nrow(indels)) {
      break; #finish
    }
    else if(indels[j,"sampleID"] != indels[(j-1),"sampleID"]){
      i <- j;
      break;
    } else {
      if(indels[j,"pos"] != (indels[(j-1),"pos"]+1)) {
        i <- j; #start a new indel
        break;
      } else {
        if(indels[j,"chr"] != indels[(j-1),"chr"]) {
          i <- j; #start a new indel
          break;
        } else {
          #This is a candidate, but check their VAFs are compatible with a Fishers exact test:
          mat <- matrix(nrow=2,ncol=2,0);
          mat[1,] <- c(indels[j-1,  "xfw"]+indels[j-1,  "xbw"], indels[j-1,  "nfw"]+indels[j-1,  "nbw"])
          mat[2,] <- c(indels[j,    "xfw"]+indels[j,    "xbw"], indels[j,    "nfw"]+indels[j,    "nbw"])
          vaf1 <- mat[1,1] / mat[1,2]
          vaf2 <- mat[2,1] / mat[2,2]
          if (abs(vaf1-vaf2) > deltavaf & mat[1,1] != mat[2,1]) {
            cat(" Breaking up indel because VAFs do not match\n");
            cat("             ",j-1, " vs ", j, ": deltavaf=", abs(vaf1-vaf2), " [",indels[j-1,"pos"],"-",indels[j,"pos"],"]",sep="");
            cat("    (mat=", mat[1,1],",",mat[1,2],",",mat[2,1],",",mat[2,2],")\n",sep="");
            i <- j;
            break;
          }
          indel_to <- indels[j,"pos"];
          to_index <- j;
          deleted_seq <- paste(deleted_seq,indels[j,"ref"],sep="");
        }
      } 
    }
  }
  #cat("   [Sample=",sample_f,"] Indel goes from=",indel_from,", to=", indel_to," [",deleted_seq,">-]\n",sep="");
  indels[c(from_index:to_index),"groupID"    ] <- group_counter;
  indels[c(from_index:to_index),"deleted_seq"] <- deleted_seq;
  group_counter <- group_counter + 1;
}
mutations$indel_group <- NA;
mutations$deleted_seq <- NA;

if(nrow(indels) > 0){
  mutations[rownames(indels),c("indel_group","deleted_seq")] = indels[rownames(indels),c("groupID","deleted_seq")]
}

mutations$qval = p.adjust(mutations$pval, method="BH", n=sites_with_100x*5)

ok_muts <- mutations[which(mutations$qval <= 0.05),]

sub <- ok_muts[grep("[-I]",ok_muts[,"mut"],invert=T),]
sub$mut <- as.character(sub$mut)

ins <- ok_muts[grep("I",   ok_muts[,"mut"]),         ]
del <- ok_muts[grep("-",   ok_muts[,"mut"]),         ]
sub <- sub[order(sub$sampleID, sub$chr, sub$pos),    ]
ins <- ins[order(ins$sampleID, ins$chr, ins$pos),    ]
del <- del[order(del$sampleID, del$chr, del$pos),    ]

new_mutations <- mutations[0,]

# Deletions (defined in mutations$indel_group):
# For every "OK" deletion, get its del-groupID and find all the other deletions belonging
# to that group. Merge them and create a new entry in mutations: combine pvalues, vaf, etc
# For every "OK" deletion first check it hasn't been already merged
for(j in 1:nrow(del)) {
  indel_group       <- del[j,"indel_group"]
  if(nrow(new_mutations[which(new_mutations$indel_group==indel_group),]) > 0) {
    next; #we already have one from the group of indels
  }
  indels_from_group                                  <- mutations[which(mutations$indel_group==indel_group),]
  new_mutations                                      <- rbind(new_mutations,indels_from_group[1,])
  new_mutations[nrow(new_mutations),"pos"          ] <- min  (indels_from_group$pos              )
  new_mutations[nrow(new_mutations),"vaf"          ] <- mean (indels_from_group$vaf              )
  new_mutations[nrow(new_mutations),"tum_globalvaf"] <- mean (indels_from_group$tum_globalvaf    )
  new_mutations[nrow(new_mutations),"pval"         ] <- min  (indels_from_group$pval             )
  new_mutations[nrow(new_mutations),"qval"         ] <- min  (indels_from_group$qval             )
  new_mutations[nrow(new_mutations),"ref"          ] <- indels_from_group[1,"deleted_seq"]
}

# Insertions. No need to look for consecutive INS. Just add them to new_mutations
new_mutations <- rbind(new_mutations,ins);

# Substitutions: merge consecutive... [using Iñigo's code]
d = sub$pos-(1:nrow(sub))
runs = rle(d)
rmpos = rep(0,nrow(sub))
runstarts = cumsum(runs$length)-runs$length+1
for (h in 1:length(runs$length)) {
  if (runs$length[h]>1) { # Adjacent mutations
    mutcluster                         = runstarts[h]:(runstarts[h]+runs$lengths[h]-1)
    if (length(unique(sub$sampleID[mutcluster]))==1 & all(abs(diff(sub$vaf[mutcluster]))<deltavaf)) {
        rmpos[mutcluster[-1]             ] = 1 # Removing all the affected rows except the first one (which we will edit to capture the complex event)
        sub[mutcluster[1],"ref"          ] = paste(sub[mutcluster,"ref"          ],collapse="")
        sub[mutcluster[1],"mut"          ] = paste(sub[mutcluster,"mut"          ],collapse="")
        sub[mutcluster[1],"mu"           ] = mean (sub[mutcluster,"mu"           ]            )
        sub[mutcluster[1],"tum_globalvaf"] = mean (sub[mutcluster,"tum_globalvaf"]            )
        sub[mutcluster[1],"vaf"          ] = mean (sub[mutcluster,"vaf"          ]            )
        sub[mutcluster[1],"pval"         ] = min  (sub[mutcluster,"pval"         ]            )
        sub[mutcluster[1],"qval"         ] = min  (sub[mutcluster,"qval"         ]            )
    }
  }
}
sub = sub[!rmpos,]        
new_mutations <- rbind(new_mutations,sub);

good_muts     <- new_mutations
good_muts$mut_site <- paste(good_muts$pos,good_muts$ref,good_muts$mut,sep = "_")
good_muts$mut_id <- paste(good_muts$sampleID,good_muts$mut_site,sep = "_")

saveRDS(object = good_muts, "shearwater_muts_all_samples_qval0.05_mincov100x_panel_1_run_1.rds")


### Code to merge the calls from both runs

# Function to intersect PASS mutations from file 1 (qval<0.05) with all mutations from file 2
relaxed_intersect = function(file1, file2, minq1 = 0.05, minp2 = 0.01) {
    
    m1 = readRDS(file1)
    m1 = m1[m1$qval<minq1, ]
    m2 = readRDS(file2)
    m2 = m2[m2$pval<minp2, ]
    
    gr1 = GenomicRanges::GRanges(m1$sampleID, IRanges::IRanges(m1$pos,m1$pos+nchar(m1$ref)-1))
    gr2 = GenomicRanges::GRanges(m2$sampleID, IRanges::IRanges(m2$pos,m2$pos+nchar(m2$ref)-1))
    ol = as.data.frame(GenomicRanges::findOverlaps(gr1, gr2, type="any", select="all"))
    m1sig = m1[unique(ol[,1]), ]
    return(m1sig)
}

run1sig = relaxed_intersect(file1 = "shearwater_muts_all_samples_qval0.05_mincov100x_panel_1_run_1.rds",
                            file2 = "initial_shearwater_calls_background_panel_1_run_2.rds")

run2sig = relaxed_intersect(file1 = "shearwater_muts_all_samples_qval0.05_mincov100x_panel_1_run_2.rds",
                            file2 = "initial_shearwater_calls_background_panel_1_run_1.rds")

# Renaming the columns of both mutation tables

run1sig$mut = as.vector(run1sig$mut)
run2sig$mut = as.vector(run2sig$mut)
run1sig$length1 = nchar(run1sig$ref)
run2sig$length2 = nchar(run2sig$ref)

run1sig = run1sig[,c("sampleID","chr","pos","ref","mut","xfw","xbw","nfw","nbw","pval","qval","length1","mut_id")]
colnames(run1sig) = c("sampleID","chr","pos","ref","mut","r1_xfw","r1_xbw","r1_nfw","r1_nbw","r1_pval","r1_qval","r1_length","mut_id")
run2sig = run2sig[,c("sampleID","chr","pos","ref","mut","xfw","xbw","nfw","nbw","pval","qval","length2","mut_id")]
colnames(run2sig) = c("sampleID","chr","pos","ref","mut","r2_xfw","r2_xbw","r2_nfw","r2_nbw","r2_pval","r2_qval","r2_length","mut_id")

# Merging both tables:
# 1. First, we merge all perfect matches.
# 2. Second, we pick the largest event when matches are imperfect.

part1 = merge(run1sig, run2sig, by=c("sampleID","chr","pos","ref","mut","mut_id")) # Merging perfect matches
part1$match = "perfect"

gr1 = GenomicRanges::GRanges(run1sig$sampleID, IRanges::IRanges(run1sig$pos,run1sig$pos+nchar(run1sig$ref)-1))
gr2 = GenomicRanges::GRanges(run2sig$sampleID, IRanges::IRanges(run2sig$pos,run2sig$pos+nchar(run2sig$ref)-1))
ol = as.data.frame(GenomicRanges::findOverlaps(gr1, gr2, type="any", select="all"))
ol1longer = ol[run1sig$r1_length[ol[,1]] > run2sig$r2_length[ol[,2]], ]
ol2longer = ol[run1sig$r1_length[ol[,1]] < run2sig$r2_length[ol[,2]], ]

# Breaking down a mut_id or vector of mut_ids into a vector of component SNVs
break_mutstr = function(str) {
    auxall = strsplit(str, split = "_")
    out = NULL # Init
    for (j in 1:length(auxall)) {
        aux = auxall[[j]]
        pos = as.numeric(aux[2])
        from = strsplit(aux[3], split="")[[1]]
        to = strsplit(aux[4], split="")[[1]]
        out = c(out, paste(seq(pos,pos+length(from)-1), from, to, sep="_")) # Original mutation breaking down into SNVs
    }
    return(out)
}

# Merging events when run1 is longer

longest_match = function(a, b, ol, aname = "run1", bname = "run2") {
    avec = unique(ol[,1]) # List of query elements
    df = NULL
    for (j in 1:length(avec)) {
        bvec = ol[ol[,1]==avec[j],2]
        astr = a$mut_id[avec[j]]
        bstr = b$mut_id[bvec]
        matchstr = paste(c(unique(paste(aname, astr, sep=":")), unique(paste(bname, bstr, sep=":"))), collapse=";")
        
        if (any(break_mutstr(bstr) %in% break_mutstr(astr))) { # We have a match
            rowj = cbind(a[avec[j],],b[bvec[1],],match=matchstr)[colnames(part1)]
            df = rbind(df,rowj)
        }
    }
    df$match = as.vector(df$match)
    return(df)
}

part2 = longest_match(a = run1sig, b = run2sig, ol = ol1longer, aname = "run1", bname = "run2")
part3 = longest_match(a = run2sig, b = run1sig, ol = ol2longer[,2:1], aname = "run2", bname = "run1")

mutations = rbind(part1, part2, part3)
mutations = mutations[order(mutations$sampleID, mutations$pos, mutations$mut), ]


