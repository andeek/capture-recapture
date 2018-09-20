# libraries + seed ----
library(foreach) # parallel for loops
library(doParallel) # parallel backend
library(LCMCR) # crc
set.seed(1234) # reproducible

# prep for parallel ----
cl <- makeCluster(10)
registerDoParallel(cl)

# load contings ----
load("results/error_simulation/contings.Rdata")

# capture-recapture ----
K <- 2^D - 1 ## from paper, number of unique capture histories

error_add_pop_N <- foreach(r = seq_along(error_add_conting), .combine = rbind) %dopar% {
  error_conting <- error_add_conting[[r]]
  `%dopar%` <- foreach::`%dopar%`
  foreach::foreach(i = seq_along(error_conting), .combine = rbind) %dopar% {
    cat(paste("error: ", error_levels[r], ", iter:", i, "\r"))
    
    conting <- error_conting[[i]][-1,]
    conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
    crc_sampler <- LCMCR::lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
    pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
    data.frame(pop_N = pop_N, i = i, r = error_levels[r])
  }
}
error_remove_pop_N <- foreach(r = seq_along(error_remove_conting), .combine = rbind) %dopar% {
  error_conting <- error_remove_conting[[r]]
  `%dopar%` <- foreach::`%dopar%`
  foreach::foreach(i = seq_along(error_conting), .combine = rbind) %dopar% {
    cat(paste("error: ", error_levels[r], ", iter:", i, "\r"))
    
    conting <- error_conting[[i]][-1,]
    conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
    crc_sampler <- LCMCR::lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
    pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
    data.frame(pop_N = pop_N, i = i, r = error_levels[r])
  }
}
# stop parallel ----
stopCluster(cl)
stopImplicitCluster()

# save ----
save(error_add_pop_N, error_remove_pop_N, file = "results/error_simulation/crc.Rdata")