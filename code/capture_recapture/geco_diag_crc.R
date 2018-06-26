## command line args ---- 
## duplication levels f; string distortion levels 5, 10, 15
# Rscript jasa_eber.R 5 5
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in the dup level (5) and dist level (5, 10, 15).", call.=FALSE)
dup_level <- as.numeric(args[1])
dist_level <- as.numeric(args[2])


## libraries ----
library(foreach) # parallel computation
library(doParallel)

get_dsn <- function(clusters) {
  ## get inclusion tables based on cluster
  db_clust <- data.frame(db = unlist(lapply(seq_along(noisy_dup_db), function(i) rep(i, nrow(noisy_dup_db[[i]])))),
                         clust = clusters)
  
  included <- data.frame(record = unique(db_clust$clust))
  for(i in seq_along(noisy_dup_db)) {
    included[, paste0("db", i)] <- as.factor(as.numeric(included$record %in% db_clust[db_clust$db == i, "clust"]))
  }
  
  as.data.frame(table(included[, -1]))
}

get_diag <- function(clusters, identity_freq) {
  `%>%` <- dplyr::`%>%`
  
  get_dsn(clusters) %>%
    dplyr::rename(Freq2 = Freq) %>%
    dplyr::left_join(identity_freq) %>%
    dplyr::mutate(diff = Freq - Freq2) %>%
    dplyr::mutate_at(dplyr::vars(dplyr::starts_with("db")), function(x) as.numeric(as.character(x))) %>%
    dplyr::mutate(num_include = db1 + db2 + db3 + db4 + db5) %>%
    dplyr::group_by(num_include) %>%
    dplyr::summarise(truth_clust_diff = sum(diff)) %>%
    dplyr::filter(num_include > 0)
}

## data & results load ----
load(paste0("data/geco_sim/geco_", dup_level, "dup_", dist_level, "dist.Rdata"))
load(paste0("results/geco_sim/eber_", dup_level, "dup_", dist_level, "dist.Rdata"))

## construct identity vector & get dsn
identity_vec <- do.call(c, identity)
identity_freq <- get_dsn(identity_vec)

## burnin + thin based on looking at diag files by hand
m <- ncol(lambda)
burnin <- 10000
thinning <- 40
lambda <- lambda[, seq_len(m) > burnin & seq_len(m) %% thinning == 0]


## mpmms ----
mpmms_linkage <- eber::smpmms.linkage(list(lambda.chain = lambda))
mpmms_lambda <- eber::clust2memb(mpmms_linkage)

mpmms_dsn <- data.frame(num_include = 1:5) # get all possibles
mpmms_dsn <- dplyr::left_join(mpmms_dsn, get_diag(mpmms_lambda, identity_freq)) # get dsn
mpmms_dsn$truth_clust_diff[is.na(mpmms_dsn$truth_clust_diff)] <- 0 # replace na with 0

## each lambda iter
## register parallel backend
cl <- detectCores()
registerDoParallel(cl)

lambda_dsn <- foreach(i = seq_len(ncol(lambda)), .combine = rbind) %dopar% {
  cat(paste("iter:", i, "\r"))
  
  lambda_dsn <- data.frame(num_include = 1:5, iter = i) # get all possibles
  lambda_dsn <- dplyr::left_join(lambda_dsn, get_diag(lambda[, i], identity_freq)) # get dsn
  lambda_dsn$truth_clust_diff[is.na(lambda_dsn$truth_clust_diff)] <- 0 # replace na with 0
  
  lambda_dsn
}


## save results
save(mpmms_dsn, lambda_dsn, file = paste0("results/geco_sim/crc_diag_", dup_level, "dup_", dist_level, "dist.Rdata"))





