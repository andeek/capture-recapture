## command line args ---- 
## duplication levels 5; string distortion levels 1, 3, 5; folder name
# Rscript jasa_eber.R 5 5
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 3) stop("Pass in the dup level (5, 10, 30, 50) and dist level (5, 10, 15).", call.=FALSE)
if (!(args[1] %in% c(5))) stop("Pass in the duplication level (5, 10, 30, 50)", call.=FALSE)
if (!(args[2] %in% c(1, 3, 5))) stop("Pass in the distortion level (5, 10, 15)", call.=FALSE)
if (!(args[3] %in% c("geco_sim"))) stop("Pass in the data folder name", call.=FALSE)
dup_level <- as.numeric(args[1])
dist_level <- as.numeric(args[2])
folder_name <- args[3]

## libraries ----
library(eber) # record linkage
set.seed(1234)

## data load ----
load(paste0("data/", folder_name, "/geco_", dup_level, "dup_", dist_level, "dist.Rdata"))

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
save(lambda, file = paste0("results/", folder_name, "/eber_", dup_level, "dup_", dist_level, "dist.Rdata"))



