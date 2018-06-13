## command line args ---- 
## duplication levels 5, 10, 30, 50; string distortion levels 5, 10, 15; folder name
# Rscript jasa_eber.R 5 5
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Pass in the dup level (5, 10, 30, 50) and dist level (5, 10, 15).", call.=FALSE)
if (!(args[1] %in% c(5, 10, 30, 50))) stop("Pass in the duplication level (5, 10, 30, 50)", call.=FALSE)
if (!(args[2] %in% c(5, 10, 15))) stop("Pass in the distortion level (5, 10, 15)", call.=FALSE)
if (!(args[3] %in% c("jasa_sim", "jasa_sim_small"))) stop("Pass in the data folder name", call.=FALSE)
dup_level <- as.numeric(args[1])
dist_level <- as.numeric(args[2])
folder_name <- args[3]


## libraries ----
library(eber) # record linkage
library(foreach) # parallel computation
library(doParallel)
library(ggplot2) # plotting
library(tidyr) # manipulation
library(dplyr)

get_diag <- function(dup_level, dist_level, nclust = 10) {
  ## data & results load ----
  load(paste0("data/", folder_name, "/jasa_", dup_level, "dup_", dist_level, "dist.Rdata"))
  load(paste0("results/", folder_name, "/eber_", dup_level, "dup_", dist_level, "dist.Rdata"))
  
  ## construct identity vector & combine data
  data <- do.call(rbind, noisy_dup_db)
  identity_vec <- do.call(c, identity)
  identity_pair <- clust2pairs(memb2clust(identity_vec))

  ## compute metrics
  ## register parallel backend
  cl <- makePSOCKcluster(nclust)
  registerDoParallel(cl)
  
  eval <- foreach(i = seq_len(ncol(lambda)), .combine = rbind) %dopar% {
    lam_i <- eber::clust2pairs(eber::memb2clust(lambda[, i]))
    
    ## precision and recall
    conf <- eber::confusion.matrix(lam_i, identity_pair)
    
    ## singletons
    single <- sum(table(lambda[, i]) == 1)
    
    ## return in df
    data.frame(tp = conf$TP, fp = conf$FP, fn = conf$FN, 
               recall = conf$TP/(conf$TP + conf$FN), precision = conf$TP/(conf$TP + conf$FP),
               num_single = single)
  }
  
  ## stop parallel backend
  stopCluster(cl)
  
  ## precision/recall trace plot
  eval %>%
    mutate(iter = 1:n()) %>%
    gather(metric, value, precision, recall) %>%
    ggplot() +
    geom_line(aes(iter, value)) +
    facet_grid(.~metric, scales = "free_y") +
    xlab("Iteration") +
    ylab("") -> p.eval
  
  ## singletons trace plot
  eval %>%
    mutate(iter = 1:n()) %>%
    ggplot() +
    geom_line(aes(iter, num_single)) +
    xlab("Iteration") +
    ylab("Number of Singletons") -> p.singleton
  
  ## singletons acf
  single_acf <- acf(eval$num_single, plot = FALSE)
  
  ## results
  res <- list(eval = eval,
              singletons = singletons,
              single_acf = single_acf,
              p.eval = p.eval,
              p.singletons = p.singleton)
  return(res)
}

## get diagnostics for specified levels ----
diag_res <- get_diag(dup_level, dist_level)

## save results
save(diag_res, file = paste0("results/", folder_name, "/eber_diag_", dup_level, "dup_", dist_level, "dist.Rdata"))







