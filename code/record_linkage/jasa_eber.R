## command line args ---- 
## duplication levels 5, 10, 30, 50; string distortion levels 5, 10, 15
# Rscript jasa_eber.R 5 5
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2) stop("Pass in the dup level (5, 10, 30, 50) and dist level (5, 10, 15).", call.=FALSE)
if (!(args[1] %in% c(5, 10, 30, 50))) stop("Pass in the duplication level (5, 10, 30, 50)", call.=FALSE)
if (!(args[2] %in% c(5, 10, 15))) stop("Pass in the distortion level (5, 10, 15)", call.=FALSE)
dup_level <- as.numeric(args[1])
dist_level <- as.numeric(args[2])

## libraries ----
library(eber) # record linkage
set.seed(1234)

## data load ----
load(paste0("data/jasa_sim/jasa_", dup_level, "dup_", dist_level, "dist.Rdata"))

## data format and params for use with eber ----
data_all <- as.matrix(do.call(rbind, noisy_dup_db))
fields <- c("year", "month", "title", "author")
string.fields <- c("title", "author")
file.start <- head(cumsum(c(1, unlist(lapply(noisy_dup_db, nrow)))), -1)
num.iter <- 30000
prior.alpha <- 1
prior.beta <- 99
pop.size <- nrow(data_all) # this might not be a good idea
num.chains <- 3

## run eber ----
out <- eber(data_all, fields, string.fields, file.start, num.iter, 
            prior.alpha, prior.beta, 
            pop.size = pop.size, num.chains = num.chains)
out.bind <- bind.chains(out)
lambda <- out.bind$lambda.chain

## save results ----
save(lambda, file = paste0("results/jasa_sim/eber_", dup_level, "dup_", dist_level, "dist.Rdata"))



