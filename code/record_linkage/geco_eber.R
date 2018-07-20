## command line args ---- 
## within db duplication levels 5; string distortion levels 5, 10, 15; num_dist 1, 2, 3
# Rscript code/data_munging/jasa_pubs_10yr.R 5 5
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Pass in the dup level and dist level.", call.=FALSE)
dup_level <- as.numeric(args[1])
dist_level <- as.numeric(args[2])
num_dist <- as.numeric(args[3])

## libraries ----
library(eber) # record linkage
set.seed(1234)

## data load ----
load(paste0("data/geco_sim/geco_", dup_level, "dup_", dist_level, "dist_", num_dist, "num.Rdata"))

## data format and params for use with eber ----
data_all <- as.matrix(do.call(rbind, noisy_dup_db))
fields <- c("fname", "lname", "by", "bm", "bd")
string.fields <- c("fname", "lname")
file.start <- head(cumsum(c(1, unlist(lapply(noisy_dup_db, nrow)))), -1)
num.iter <- 50000
prior.alpha <- 1
prior.beta <- 99
pop.size <- nrow(data_all) # this might not be a good idea
num.chains <- 1

## run eber ----
out <- eber(data_all, fields, string.fields, file.start, num.iter, 
            prior.alpha, prior.beta, 
            pop.size = pop.size, num.chains = num.chains)
out.bind <- bind.chains(out)
lambda <- out.bind$lambda.chain

## save results ----
save(lambda, file = paste0("data/geco_sim/geco_", dup_level, "dup_", dist_level, "dist_", num_dist, "num.Rdata"))



