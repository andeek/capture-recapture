# libraries + seed ----
library(foreach) # parallel for loops
library(doParallel) # parallel backend
library(LCMCR) # crc
set.seed(1234) # reproducible

# global params ----
rep <- 100

# load true conting ----
load("results/error_simulation/contings.Rdata")

# keep singletons, but randomly choose multis ----
# copy of true_conting with labels
tab <- true_conting
tab[-ncol(tab)] <- apply(tab[-ncol(tab)], 2, function(x) as.logical(as.character(x))) # change factors to logical
tab$type <- ifelse(rowSums(tab[, -ncol(tab)]) == 0, "zero", ifelse(rowSums(tab[, -ncol(tab)]) == 1, "single", "multi")) # classify buckets

# count the number of buckets and the number of items
num_multi <- sum(tab[tab$type == "multi", "Freq"])
num_multi_buckets <- sum(tab$type == "multi")

# randomly select from multinomial with equal probs
new_multi <- rmultinom(rep, num_multi, rep(1, num_multi_buckets))
contings <- apply(new_multi, 2, function(multi) { tab1 <- tab; tab1[tab1$type == "multi", "Freq"] <- multi; tab1[, -ncol(tab1)] })

# prep for parallel ----
cl <- makeCluster(16)
registerDoParallel(cl)

# capture-recapture ----
pop_N <- foreach(i = seq_along(contings), .combine=rbind, .verbose = TRUE) %dopar% {
    conting <- error_add_conting[[i]][-1,]
    K <- sum(conting$Freq > 0) ## from paper, number of unique capture histories
    conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
    crc_sampler <- LCMCR::lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
    pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
    data.frame(pop_N = pop_N, i = i)
  }

# save ----
save(pop_N, file = "results/error_simulation/crc_multi.Rdata")


