library("GenomicRanges")
library("deepSNV")
library("Rsamtools")

args <- commandArgs(TRUE)
background_panel <- args[1]
run <- args[2]
sample_start <- args[3]
sample_end <- args[4]


sample_list <- read.table("sample_list_1181.tsv", sep = "\t", stringsAsFactors = F, header = F)

genome_file = "MN908947.3.fasta"

logbb = deepSNV:::logbb
dbetabinom = VGAM::dbetabinom

min_totest = 2

ref_list <- c("A","T","C","G","-","INS")

normcounts_obj <- list()
normcounts_obj$counts <- readRDS(file = paste0("Analysis/Shearwater/background_panel_",background_panel,".rds"))
norm_total <- apply(normcounts_obj$counts[, , 1:6] + normcounts_obj$counts[, , 7:12], c(2, 3), sum)
norm_total <- norm_total / rowSums(norm_total)

normal_rho_panel <- readRDS(paste0("Analysis/Shearwater/background_panel_",background_panel,"_rho_mu_calculation.rds"))
normal_rho_panel$rho[which(normal_rho_panel$rho <= 1e-03)] <- 1e-03

tumcounts_obj_all <- list()
tumcounts_obj_all$counts <- readRDS(file = paste0("Data/run_",run,"_coverage.rds"))

tum_total <- apply(tumcounts_obj_all$counts[, , 1:6] + tumcounts_obj_all$counts[, , 7:12], c(2, 3), sum)
tum_total <- tum_total / rowSums(tum_total)

mutations_allsamples = NULL

for (k in sample_start:sample_end) {
  tumcounts_obj = tumcounts_obj_all
  tumcounts_obj$counts = array(tumcounts_obj$counts[k, , ], dim = c(1, dim(tumcounts_obj$counts[k, , ])))
  
  mutsites = which((tumcounts_obj$counts[1, , 1:6] + tumcounts_obj$counts[1, , 7:12]) >= min_totest & norm_total < 0.95, arr.ind = T)
  
  if (nrow(mutsites) > 0) {
    mutations = data.frame(
      sampleID = sample_list[k,],
      chr = "MN908947.3",
      pos = c(1:29903)[mutsites[, 1]],
      ref = NA,
      mut = ref_list[mutsites[, 2]],
      xfw = NA,
      xbw = NA,
      nfw = NA,
      nbw = NA,
      mu = NA,
      pval = NA,
      rho = NA
    )
    mutations$tum_globalvaf = apply(mutsites, 1, function(x) tum_total[x[1], x[2]])
    
    l = 6
    indsm = cbind(mutsites[, 2], mutsites[, 2] + l)
    sites_gr = GRanges(mutations$chr, IRanges(mutations$pos, mutations$pos))
    seqs = as.character(scanFa(genome_file, sites_gr))
    mutations$ref = as.character(seqs)
    
    for(x in 1:nrow(mutations)){
      inds = indsm[x, ]
      norcounts = normcounts_obj$counts[, mutsites[x], ]
      tumcounts = tumcounts_obj$counts[, mutsites[x], ]
      
      pseudo = .Machine$double.eps
      x.fw = tumcounts[inds[1]]
      x.bw = tumcounts[inds[2]]
      x.tot = x.fw + x.bw
      
      n.fw = sum(tumcounts[1:l])
      n.bw = sum(tumcounts[(l + 1):(2 * l)])
      n.tot = n.fw + n.bw
      
      X.fw = sum(norcounts[, inds[1]])
      X.bw = sum(norcounts[, inds[2]])
      X.tot = X.fw + X.bw
      
      N.fw = sum(norcounts[, 1:l])
      N.bw = sum(norcounts[, (l + 1):(2 * l)])
      N.tot = N.fw + N.bw
      
      mu = max(X.tot, pseudo) / max(N.tot, pseudo)
      
      counts = cbind(rowSums(norcounts[, 1:l]), rowSums(norcounts[, (l +
                                                                       1):(2 * l)]), norcounts[, inds])
      rho = normal_rho_panel$rho[which(normal_rho_panel$pos == mutsites[x,1] & normal_rho_panel$mut == ref_list[mutsites[x,2]])]
      
      rdisp = (1 - rho) / rho
      
      prob0.tot = (X.tot + x.tot) / (N.tot + n.tot)
      prob0.tot[prob0.tot == 0] = pseudo
      prob1s.tot = x.tot / (n.tot + pseudo)
      prob1s.tot[prob1s.tot == 0] = pseudo
      prob1c.tot = X.tot / (N.tot + pseudo)
      prob1c.tot[prob1c.tot == 0] = pseudo
      prob1s.tot = pmax(prob1s.tot, prob1c.tot) # Min error rate is that of the population (one-sided test)
      # Modified by fa8 to avoid p=1, which results in beta=0
      nu0.tot = prob0.tot * rdisp
      nu1s.tot = min(1 - pseudo, prob1s.tot) * rdisp
      nu1c.tot = min(1 - pseudo, prob1c.tot) * rdisp
      
      LL.tot = logbb(x.tot, n.tot, nu0.tot, rdisp) + logbb(X.tot, N.tot, nu0.tot, rdisp) - logbb(x.tot, n.tot, nu1s.tot, rdisp) - logbb(X.tot, N.tot, nu1c.tot, rdisp)
      pvals_tot = pchisq(-2 * LL.tot, df = 1, lower.tail = F) / 2 # We divide by 2 as we are performing a 1-sided test
      mutations[x, 6:12] = c(x.fw, x.bw, n.fw, n.bw, mu, pvals_tot, rho)
    }
    mutations_allsamples = rbind(mutations_allsamples, mutations)
  }
}

saveRDS(object = mutations_allsamples, file = paste0("Analysis/Shearwater/background_panel_",background_panel,"_run_",run,"/initial_shearwater_calls_background_panel_",background_panel,"_run_",run,"_samples_",sample_start,"_",sample_end,".rds"))
