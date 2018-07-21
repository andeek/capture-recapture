## command line args ---- 
## duplication levels 5; string distortion levels 5, 10, 15;
# Rscript jasa_eber.R 5 5
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Pass in the dup level (5) and dist level (5, 10, 15).", call.=FALSE)
dup_level <- as.numeric(args[1])
dist_level <- as.numeric(args[2])
num_dist <- as.numeric(args[3])


## libraries ----
library(dga) # crc

set.seed(42)

## data + results load ----
load(paste0("data/geco_sim/geco_", dup_level, "dup_", dist_level, "dist_", num_dist, "num.Rdata"))
load(paste0("results/geco_sim/eber_", dup_level, "dup_", dist_level, "dist_", num_dist, "num.Rdata"))

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
  # load the graphs to make the estimates
  data(graphs5)
  
  # select expansion factor defining the largest number of unrecorded elements.
  # this makes Nmissing <- 0:(sum(Y)*fac)
  fac <- 5
  
  # set prior
  delta <- 1/2^length(noisy_dup_db)
  
  # average over all decomposible graphical models for p lists
  Nmissing <- 0:(nrow(included)*fac)
  weights <- bma.cr(table(included[, -1]), Nmissing, delta, graphs5)
  
  sample(nrow(included) + Nmissing, 10000, prob = colSums(weights), replace = TRUE)
}

## truth ----
pop_N_truth <- get_pop_N(do.call(c, identity))

## mpmms ----
mpmms_linkage <- eber::smpmms.linkage(list(lambda.chain = lambda))
mpmms_lambda <- eber::clust2memb(mpmms_linkage)
pop_N_mpmms <- get_pop_N(mpmms_lambda)

## fully bayes ----
pop_N_bayes <- apply(lambda, 2, get_pop_N)

save(pop_N_truth, mpmms_lambda, pop_N_mpmms, pop_N_bayes, file = paste0("results/geco_sim/bma_", dup_level, "dup_", dist_level, "dist_", num_dist, "num.Rdata"))




