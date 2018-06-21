## command line args ---- 
## duplication levels 5; string distortion levels 1, 5, 10;
# Rscript jasa_eber.R 5 5
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in the dup level (5) and dist level (1, 5, 10).", call.=FALSE)
dup_level <- as.numeric(args[1])
dist_level <- as.numeric(args[2])


## libraries ----
library(LCMCR) # crc

set.seed(42)

## data + results load ----
load(paste0("data/geco_sim/geco_", dup_level, "dup_", dist_level, "dist.Rdata"))
load(paste0("results/geco_sim/eber_", dup_level, "dup_", dist_level, "dist.Rdata"))

## burnin + thin based on looking at diag files by hand
m <- ncol(lambda)
burnin <- 10000
thinning <- 40
lambda <- lambda[, seq_len(m) > burnin & seq_len(m) %% thinning == 0]

## perform CRC for various "clusterings" ----
get_pop_N <- function(clusters) {
  ## get inclusion tables based on cluster
  db_clust <- data.frame(db = unlist(lapply(seq_along(noisy_dup_db), function(i) rep(i, nrow(noisy_dup_db[[i]])))),
                         clust = clusters)
  
  included <- data.frame(record = unique(db_clust$clust))
  for(i in seq_along(noisy_dup_db)) {
    included[, paste0("db", i)] <- as.factor(as.numeric(included$record %in% db_clust[db_clust$db == i, "clust"]))
  }
  
  ## capture-recapture
  K <- nrow(unique(included[, -1])) ## from paper, number of unique capture histories
  crc_sampler <- lcmCR(included[, -1], K = K, seed = 42)
  pop_N <- lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 10000, thinning = 20)
  
  pop_N
}

## truth ----
pop_N_truth <- get_pop_N(do.call(c, identity))

## mpmms ----
mpmms_linkage <- eber::smpmms.linkage(list(lambda.chain = lambda))
mpmms_lambda <- eber::clust2memb(mpmms_linkage)
pop_N_mpmms <- get_pop_N(mpmms_lambda)

## fully bayes ----
pop_N_bayes <- apply(lambda, 2, get_pop_N)

save(pop_N_truth, mpmms_lambda, pop_N_mpmms, pop_N_bayes, file = paste0("results/geco_sim/crc_", dup_level, "dup_", dist_level, "dist.Rdata"))




