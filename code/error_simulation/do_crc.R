# libraries + seed ----
library(foreach) # parallel for loops
library(doParallel) # parallel backend
library(LCMCR) # crc
set.seed(1234) # reproducible

# prep for parallel ----
cl <- makeCluster(3)
registerDoParallel(cl)

# load contings ----
load("results/error_simulation/contings.Rdata")

# capture-recapture ----
K <- 2^D - 1 ## from paper, number of unique capture histories

error_add_pop_N <- foreach(r = error_levels, .combine = rbind) %dopar% {
  error_conting <- error_add_conting[[r]]
  foreach(i = seq_along(error_conting), .combine = rbind) %dopar% {
    conting <- error_conting[[i]][-1,]
    conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
    crc_sampler <- lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
    pop_N <- lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
    data.frame(pop_N = pop_N, i = i, r = r)
  }
}
error_remove_pop_N <- foreach(r = error_levels, .combine = rbind) %dopar% {
  error_conting <- error_remove_conting[[r]]
  foreach(i = seq_along(error_conting), .combine = rbind) %dopar% {
    conting <- error_conting[[i]][-1,]
    conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
    crc_sampler <- lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
    pop_N <- lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
    data.frame(pop_N = pop_N, i = i, r = r)
  }
}

# stop parallel ----
stopCluster(cl)
stopImplicitCluster()

# save ----
save(error_add_pop_N, error_remove_pop_N, file = "results/error_simulation/crc.Rdata")