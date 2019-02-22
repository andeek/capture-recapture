## command line args ---- 
## duplication levels 5; string distortion levels 5, 10, 15;
# Rscript jasa_eber.R 5 5
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in the error_type ('remove' or 'add') and bucket_type ('single' or 'multi')", call.=FALSE)
error_type <- args[1]
bucket_type <- args[2]

# libraries + seed ----
library(foreach) # parallel for loops
library(doParallel) # parallel backend
library(LCMCR) # latent class crc
library(Rcapture) # frequentist crc
library(conting) # King BMA crc
library(dga) # BMA crc
set.seed(1234) # reproducible

# load true conting + global params ----
load("results/error_simulation/contings.Rdata")

# distort ----
# function for distorting a contingency table 
# adds or removes a fixed number of singletons or multis
# input: conting - a contingency table (df), r - a proportion to be distorted, error_type - "remove" or "add", bucket_type - "single" or "multi"
# output: distorted contingency table
distort_fixed_conting <- function(conting, r, error_type = "remove", bucket_type = "single") {
  stopifnot(error_type %in% c("remove", "add"))
  stopifnot(is.numeric(r) & r <= 1 & r >= 0)
  stopifnot(bucket_type %in% c("single", "multi"))
  
  tab <- conting # make copy
  tab[-ncol(tab)] <- apply(tab[-ncol(tab)], 2, function(x) as.logical(as.character(x))) # change factors to logical
  tab$type <- ifelse(rowSums(tab[, -ncol(tab)]) == 0, "zero", ifelse(rowSums(tab[, -ncol(tab)]) == 1, "single", "multi")) # classify buckets
  
  num_error <- ceiling(sum(tab[tab$type == "single", "Freq"])*r) # number of records to remove or add
  
  # choose where to take or add the records from
  where_type <- table(sample(which(tab$type == bucket_type), num_error, replace = TRUE))
  
  if(error_type == "remove") {
    tab[names(where_type), "Freq"] <- tab[names(where_type), "Freq"] - where_type
    
    # fix less than 0
    # redistribute
    while(sum(tab$Freq < 0) > 0) {
      where_type <- table(sample(which(tab$type == "single" & tab$Freq > 0), sum(-tab[tab$Freq < 0, "Freq"]), replace = TRUE))
      tab[tab$Freq < 0, "Freq"] <- 0
      tab[names(where_type), "Freq"] <- tab[names(where_type), "Freq"] - where_type
    }
  } else if(error_type == "add") {
    tab[names(where_type), "Freq"] <- tab[names(where_type), "Freq"] + where_type
  }
  return(tab[, -ncol(tab)])  
}

conting_tbl <- lapply(error_levels, function(r) lapply(seq_len(rep), function(i) distort_fixed_conting(true_conting, r, error_type, bucket_type)))
names(conting_tbl) <- error_levels

# prep for parallel ----
cl <- makeCluster(16)
registerDoParallel(cl)

# function for crc
do_all_crc <- function(conting) {
  # manipulations of conting
  num_conting <- conting # numeric version of conting instead of T/F
  num_conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
  non_df_conting <- as.data.frame(do.call(rbind, apply(conting, 1, function(x) matrix(rep(x[-length(x)], x[length(x)]), ncol=length(x) - 1, nrow=x[length(x)], byrow = TRUE))))
  non_df_conting <- table(non_df_conting)
  na_conting <- rbind(conting, data.frame(db1 = FALSE, db2 = FALSE, db3 = FALSE, db4 = FALSE, db5 = FALSE, Freq = NA))
  
  # latent class crc
  K <- sum(conting$Freq > 0) ## from paper, number of unique capture histories
  lc_crc_sampler <- LCMCR::lcmCR(num_conting, tabular = TRUE, K = K, seed = 1234)
  lc_pop_N <- LCMCR::lcmCR_PostSampl(lc_crc_sampler, burnin = 100000, samples = 500, thinning = 20)
  lc_pop_N_ci <- quantile(lc_pop_N, probs = c(.025, .975))
  
  # BMA crc
  data(graphs5) # load the graphs to make the estimates
  fac <- 5 # select expansion factor defining the largest number of unrecorded elements. this makes Nmissing <- 0:(sum(Y)*fac)
  n <- sum(conting[, ncol(conting)])
  delta <- 1/2^(ncol(conting) - 1) # set prior
  Nmissing <- 0:(n*fac) # average over all decomposible graphical models for p lists
  weights <- bma.cr(non_df_conting, Nmissing, delta, graphs5)
  bma_pop_N <- sample(n + Nmissing, 10000, prob = colSums(weights), replace = TRUE)
  bma_pop_N_ci <- quantile(bma_pop_N, probs = c(.025, .975))
  
  # King BMA crc
  mcmc_iter <- 55000
  king_crc_model <- bict(formula = Freq ~ (db1 + db2 + db3 + db4 + db5)^4, data = na_conting, n.sample = mcmc_iter)
  N <- total_pop(king_crc_model)$TOT
  burnin <- 5000 # get posterior results
  thin <- 50
  idx <- seq_len(mcmc_iter) > burnin & seq_len(mcmc_iter) %% thin == 0
  king_pop_N <- N[idx]
  king_pop_N_ci <- quantile(king_pop_N, probs = c(.025, .975))
  
  # frequentist crc
  freq_pop_N <- closedpCI.t(apply(conting, 2, function(x) as.integer(x)), dfreq = TRUE, m = "Mth")
  freq_pop_N_ci <- freq_pop_N$CI[, c("infCL", "supCL")]
  
  # results
  data.frame(matrix(c(lc_pop_N_ci, bma_pop_N_ci, king_pop_N_ci, freq_pop_N_ci), ncol = 2, byrow = TRUE), 
             method = c("lc", "bma", "king", "freq")) -> res
  names(res)[1:2] <- c("lb_95", "ub_95")
  
  return(res)
}

pop_N_ci <- foreach(r = seq_along(error_levels), .combine=rbind, .verbose = TRUE) %:% 
  foreach(i = seq_len(rep), .combine=rbind, .verbose = TRUE) %dopar% {
    conting <- conting_tbl[[r]][[i]][-1,]
    res <- do_all_crc(conting)
    data.frame(res, true_N = sum(single_remove_conting[[r]][[i]][,"Freq"]), i = i, r = error_levels[r])
  }

save(pop_N_ci, file = paste0("results/error_simulation/crc_fixed_multi_", error_type, "_", bucket_type, ".Rdata"))
