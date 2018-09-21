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

# permute multis
multi_conting <- lapply(seq_len(rep), function(i) {
  tab <- true_conting
  tab[-ncol(tab)] <- apply(tab[-ncol(tab)], 2, function(x) as.logical(as.character(x))) # change factors to logical
  tab$type <- ifelse(rowSums(tab[, -ncol(tab)]) == 0, "zero", ifelse(rowSums(tab[, -ncol(tab)]) == 1, "single", "multi")) # classify buckets
  
  # permute multis
  tab[which(tab$type == "multi"), "Freq"] <- tab[sample(which(tab$type == "multi"), sum(tab$type == "multi")), "Freq"]
  tab[, -ncol(tab)]
})
all_conting <- lapply(seq_len(rep), function(i) {
  tab <- true_conting
  tab[-ncol(tab)] <- apply(tab[-ncol(tab)], 2, function(x) as.logical(as.character(x))) # change factors to logical
  tab$type <- ifelse(rowSums(tab[, -ncol(tab)]) == 0, "zero", ifelse(rowSums(tab[, -ncol(tab)]) == 1, "single", "multi")) # classify buckets
  
  # permute multis
  tab[which(tab$type != "zero"), "Freq"] <- tab[sample(which(tab$type != "zero"), sum(tab$type != "zero")), "Freq"]
  tab[, -ncol(tab)]
})

# prep for parallel ----
cl <- makeCluster(16)
registerDoParallel(cl)

# capture-recapture ----
multi_pop_N <- foreach(i = seq_along(multi_conting), .combine=rbind, .verbose = TRUE) %dopar% {
    conting <- multi_conting[[i]][-1,]
    K <- sum(conting$Freq > 0) ## from paper, number of unique capture histories
    conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
    crc_sampler <- LCMCR::lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
    pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
    data.frame(pop_N = pop_N, i = i)
}
all_pop_N <- foreach(i = seq_along(all_conting), .combine=rbind, .verbose = TRUE) %dopar% {
  conting <- all_conting[[i]][-1,]
  K <- sum(conting$Freq > 0) ## from paper, number of unique capture histories
  conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
  crc_sampler <- LCMCR::lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
  pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
  data.frame(pop_N = pop_N, i = i)
}

pop_N <- rbind(data.frame(multi_pop_N, type = "multi"), data.frame(all_pop_N, type = "all"))

# save ----
save(pop_N, file = "results/error_simulation/crc_multi.Rdata")


