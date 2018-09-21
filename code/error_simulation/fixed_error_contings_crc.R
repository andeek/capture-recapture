# libraries + seed ----
library(foreach) # parallel for loops
library(doParallel) # parallel backend
library(LCMCR) # crc
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

single_remove_conting <- lapply(error_levels, function(r) lapply(seq_len(rep), function(i) distort_fixed_conting(true_conting, r, "remove", "single")))
single_add_conting <- lapply(error_levels, function(r) lapply(seq_len(rep), function(i) distort_fixed_conting(true_conting, r, "add", "single")))
multi_remove_conting <- lapply(error_levels, function(r) lapply(seq_len(rep), function(i) distort_fixed_conting(true_conting, r, "remove", "multi")))
multi_add_conting <- lapply(error_levels, function(r) lapply(seq_len(rep), function(i) distort_fixed_conting(true_conting, r, "add", "multi")))
names(single_remove_conting) <- names(single_add_conting) <- names(multi_remove_conting) <- names(multi_add_conting) <- error_levels

# prep for parallel ----
cl <- makeCluster(16)
registerDoParallel(cl)

# capture-recapture ----
single_remove_pop_N <- foreach(r = seq_along(error_levels), .combine=rbind, .verbose = TRUE) %:% 
  foreach(i = seq_len(rep), .combine=rbind, .verbose = TRUE) %dopar% {
    conting <- single_remove_conting[[r]][[i]][-1,]
    K <- sum(conting$Freq > 0) ## from paper, number of unique capture histories
    conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
    crc_sampler <- LCMCR::lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
    pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
    data.frame(pop_N = pop_N, i = i, r = error_levels[r])
  }
single_add_pop_N <- foreach(r = seq_along(error_levels), .combine=rbind, .verbose = TRUE) %:% 
  foreach(i = seq_len(rep), .combine=rbind, .verbose = TRUE) %dopar% {
    conting <- single_add_conting[[r]][[i]][-1,]
    K <- sum(conting$Freq > 0) ## from paper, number of unique capture histories
    conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
    crc_sampler <- LCMCR::lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
    pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
    data.frame(pop_N = pop_N, i = i, r = error_levels[r])
  }
multi_remove_pop_N <- foreach(r = seq_along(error_levels), .combine=rbind, .verbose = TRUE) %:% 
  foreach(i = seq_len(rep), .combine=rbind, .verbose = TRUE) %dopar% {
    conting <- multi_remove_conting[[r]][[i]][-1,]
    K <- sum(conting$Freq > 0) ## from paper, number of unique capture histories
    conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
    crc_sampler <- LCMCR::lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
    pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
    data.frame(pop_N = pop_N, i = i, r = error_levels[r])
  }
multi_add_pop_N <- foreach(r = seq_along(error_levels), .combine=rbind, .verbose = TRUE) %:% 
  foreach(i = seq_len(rep), .combine=rbind, .verbose = TRUE) %dopar% {
    conting <- multi_add_conting[[r]][[i]][-1,]
    K <- sum(conting$Freq > 0) ## from paper, number of unique capture histories
    conting[, -ncol(conting)] <- as.data.frame(lapply(conting[, -ncol(conting)], function(x) as.factor(as.character(as.numeric(x)))), stringsAsFactors = TRUE)
    crc_sampler <- LCMCR::lcmCR(conting, tabular = TRUE, K = K, seed = 1234)
    pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
    data.frame(pop_N = pop_N, i = i, r = error_levels[r])
  }
pop_N <- rbind(data.frame(single_add_pop_N, type = "add", bucket_type = "single"), 
               data.frame(single_remove_pop_N, type = "remove", bucket_type = "single"),
               data.frame(multi_add_pop_N, type = "add", bucket_type = "single"), 
               data.frame(multi_remove_pop_N, type = "remove", bucket_type = "single"))

# stop parallel ----
stopCluster(cl)
stopImplicitCluster()

# truth
K <- sum(true_conting$Freq > 0)
crc_sampler <- LCMCR::lcmCR(true_conting[-1,], tabular = TRUE, in_list_label = "TRUE", not_in_list_label = "FALSE", K = K, seed = 1234)
true_pop_N <- LCMCR::lcmCR_PostSampl(crc_sampler, burnin = 100000, samples = 500, thinning = 20)
pop_N <- rbind(data.frame(error_add_pop_N, type = "add"), data.frame(error_remove_pop_N, type = "remove"))

# save ----
save(pop_N, true_pop_N, file = "results/error_simulation/fixed_crc.Rdata")


