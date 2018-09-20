# libraries + seed ----
library(foreach) # parallel for loops
library(doParallel) # parallel backend
library(LCMCR) # crc
set.seed(1234) # reproducible

# prep for parallel ----
cl <- makeCluster(16)
registerDoParallel(cl)

# load contings ----
load("results/error_simulation/contings.Rdata")

# capture-recapture ----
error_add_pop_N <- foreach(r = seq_along(error_levels), .combine=rbind, .verbose = TRUE) %:% 
  foreach(i = seq_len(rep), .combine=rbind, .verbose = TRUE) %dopar% {
    conting <- error_add_conting[[r]][[i]][-1,]
    K <- sum(conting$Freq > 0) ## from paper, number of unique capture histories
    conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
    crc_sampler <- LCMCR::lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
    pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
    data.frame(pop_N = pop_N, i = i, r = error_levels[r])
  }

error_remove_pop_N <- foreach(r = seq_along(error_levels), .combine=rbind, .verbose = TRUE) %:% 
  foreach(i = seq_len(rep), .combine=rbind, .verbose = TRUE) %dopar% {
    conting <- error_remove_conting[[r]][[i]][-1,]
    K <- sum(conting$Freq > 0) ## from paper, number of unique capture histories
    conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
    crc_sampler <- LCMCR::lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
    pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
    data.frame(pop_N = pop_N, i = i, r = error_levels[r])
  }

# stop parallel ----
stopCluster(cl)
stopImplicitCluster()

# truth
K <- sum(true_conting$Freq > 0)
crc_sampler <- LCMCR::lcmCR(true_conting[-1,], tabular = TRUE, in_list_label = "TRUE", not_in_list_label = "FALSE", K = K, seed = 1234)
true_pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)

# save ----
save(error_add_pop_N, error_remove_pop_N, true_pop_N, file = "results/error_simulation/crc.Rdata")
