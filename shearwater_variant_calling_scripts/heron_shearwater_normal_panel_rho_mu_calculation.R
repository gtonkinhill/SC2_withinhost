library("GenomicRanges")
library("deepSNV")
library("Rsamtools")

logbb = deepSNV:::logbb
dbetabinom = VGAM::dbetabinom

prior        <- NULL

estimateRho_gridml = function(x, mu, prior = NULL) {
  # By Tim Coorens
  rhovec = 10 ^ seq(-6, -0.5, by = 0.05) # rho will be bounded within 1e-6 and 0.32
  mm = c(x[, 3:4])
  cov = c(x[, 1:2])
  if (is.null(prior)) {
    ll = sapply(rhovec, function(rhoj)
      sum(dbetabinom(
        x = mm,
        size = cov,
        rho = rhoj,
        prob = mu,
        log = T
      )))
  }
  else{
    ll = sapply(rhovec, function(rhoj)
      sum(dbetabinom(
        x = mm,
        size = cov,
        rho = rhoj,
        prob = mu,
        log = T
      ))) + log(prior)
  }
  rhovec[ll == max(ll)][1]
}

setwd("/Users/al28/Desktop/ProjectHeron/")

genome_file = "/Users/al28/Desktop/ProjectHeron/MN908947.3.fasta"

ref_list <- c("A","T","C","G","-","INS")

normcounts_obj <- list()
normcounts_obj$counts <- readRDS(file = "background_panel_2.rds")
norm_total <- apply(normcounts_obj$counts[, , 1:6] + normcounts_obj$counts[, , 7:12], c(2, 3), sum)
norm_total <- norm_total / rowSums(norm_total)

mutsites = array(data = NA, dim = c(29903 * 6,2))
mutsites[,1] <- rep(1:29903, times = 6)
mutsites[,2] <- rep(1:6, each = 29903)

mutations = data.frame(
  chr = "MN908947.3",
  pos = c(1:29903)[mutsites[, 1]],
  ref = NA,
  mut = ref_list[mutsites[, 2]],
  mu = NA,
  rho = NA)

l = 6
indsm = cbind(mutsites[, 2], mutsites[, 2] + l)
sites_gr = GRanges(mutations$chr, IRanges(mutations$pos, mutations$pos))
seqs = as.character(scanFa(genome_file, sites_gr))
mutations$ref = as.character(seqs)

for(x in 1:nrow(mutations)){
  inds = indsm[x, ]
  norcounts = normcounts_obj$counts[, mutsites[x], ]
  
  pseudo = .Machine$double.eps
  
  X.fw = sum(norcounts[, inds[1]])
  X.bw = sum(norcounts[, inds[2]])
  N.fw = sum(norcounts[, 1:l])
  N.bw = sum(norcounts[, (l + 1):(2 * l)])
  mu = max(X.fw + X.bw, pseudo) / max(N.fw + N.bw, pseudo)
  
  counts = cbind(rowSums(norcounts[, 1:l]), rowSums(norcounts[, (l +
                                                                   1):(2 * l)]), norcounts[, inds])
  rho = estimateRho_gridml(counts, mu, prior)
  
  mutations[x, 5:6] = c(mu, rho)
}
saveRDS(mutations, "background_panel_2_rho_mu_calculation.rds")
