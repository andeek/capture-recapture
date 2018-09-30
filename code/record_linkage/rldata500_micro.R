## libraries ----
library(MASS)
library(gtools)
library(partitionsSPscale) # record linkage model
library(RecordLinkage)
set.seed(1234)

## data load ----
data("RLdata500")

## source functions ----
source("code/record_linkage/partitions_slice_parchap.R")

## data format and params for use with partitionsSPscale ----
x <- RLdata500[, -c(2, 4)]
id <- identity.RLdata500
N <- nrow(x)
cat_fields <- c(3, 4, 5)
str_fields <- c(1, 2)
n_fields <- length(cat_fields) + length(str_fields)

# for each cat field, remap values to a consecutive list of integers starting with 1
for(cat in cat_fields) {
  x[, cat] <- as.numeric(as.factor(x[, cat]))
}

# create initial cluster assignment according
z_init <- rep(0, N)
z_init <- sample(seq_len(N), N, rep=TRUE)

# no gaps in cluster assigment #
z_init <- as.numeric(factor(z_init, labels = seq_len(length(unique(z_init)))))
z <- z_init

# params
lam_prior <- "NBD"
nsamples <- 10
burn <- 10 # burn-in period
spacing <- 10 # thin for Web_Sampler
thin1 <- thin <- 1

# beta prior for distortions
abeta <- function(bmean, bsd){
  return (bmean*(1 - bmean)*(bmean/bsd^2) - bmean)
}

# hyperparams
bmean <- 0.05 # mean
bsd <- 0.01   # sd
cb <- abeta(bmean, bsd)
db <- cb/bmean - cb
pdistortion <- cb/(cb + db)
betas <- rep(pdistortion, n_fields)

# get empirical distribution
calcProp <- function(i, data) {
  tab <- table(data[,i])
  return(tab/sum(tab))
}
proportions <- lapply(seq_len(n_fields), calcProp, data=x)

# precalculate string distances by column
calc_string_distances <- function(l, alpha, d, c) {
  alpha_l <- matrix(alpha[[l]], ncol = 1)
  S <- names(alpha[[l]])
  string_distances <- exp(-c * d(S, S))
  numerator <- matrix(rep(alpha_l, length(alpha_l)), nrow = length(alpha_l), byrow = TRUE)*string_distances
  normalized <- numerator/matrix(rep(colSums(numerator), length(alpha_l)), nrow = length(alpha_l), byrow = TRUE)
  
  rownames(normalized) <- S
  colnames(normalized) <- S
  return(normalized)
}
d <- function(string1, string2){ adist(string1, string2) }
str_proportions <- lapply(str_fields, calc_string_distances, alpha = proportions, d = d, c = 1)
# each column of str_proportions[[l]] is a distribution str_proportions[[l]][w, y]

M <- 100 # trucation for mu parameters in NBD model

# truncate the interval for slice sampling for the distortion probabilities
lods <- 0 # lower bound on support of the distribution of distortion parameters
upds <- 5/100 # upper bound on support of the distribution of distortion parameters

# labeling
rep <- "RLdata500"

# prepare for model
if(lam_prior == "PY") {
  # pitman-yor prior
  samind <- c(1, 1, 1) ## indicator vector: which parameters to slice sample
  
  theta <- N/2 # 100#concentration parameter
  delta <- 0.5 # discount parameter
  param <- c(theta, delta)
  
  # Gamma(ag, bg) prior over param theta
  bg1 <- 1/(N/2) #1/100
  ag1 <- 1
  
  # Beta(ab,bb)prior over delta
  ab1 <- 1
  bb1 <- 1
  hpriorpar <- c(ag1, bg1, ab1, bb1)
  
  # Beta prior for distortion parameters (betas)
  hpriords <- c(cb, db)
  
  w <- 1  # Size of the steps for creating interval for Slice sampler (default 1)
  m <- 10 # Limit on steps for Slice sampler(default infinite)
  lo <- c(0, 0) # Lower bound on support of the distribution of prior parameters
  up <- c(Inf, 1) # Upper bound on support of the distribution of prior parameters
  x1 <- param
  
  init <- c(theta, delta, betas)
} else if(lam_prior == "DP") {
  # dirichlet process
  samind <- c(1, 1) ## indicator vector: which parameters to slice sample
  
  theta <- N/2 #concentration parameter
  param <- as.vector(theta)
  
  # Gamma(ag1, bg1) prior over param theta
  ag1 <- 1
  bg1 <- 1/(N/2) #1/100
  hpriorpar <- c(ag1, bg1)
  
  # Beta prior for distortion parameters (betas)
  hpriords <- c(cb, db)
  
  w <- 1 #3,10 # Size of the steps for creating interval for Slice sampler (default 1)
  m <- 10 #30,100 # Limit on steps for Slice sampler(default infinite)
  lo <- 0 # Lower bound on support of the distribution of prior parameters
  up <- Inf # Upper bound on support of the distribution of prior parameters
  x1 <- param
  init <- c(x1, betas)
} else if(lam_prior == "NBNBF") {
  # NBNB microclust prior
  samind <- c(0, 0, 1, 1, 1) ## indicator vector: which parameters to slice sample
  
  anb <- 1
  qnb <- 1-2/N
  rnb <- 5
  pnb <- 0.15
  param <- c(anb, qnb, rnb, pnb)
  
  #Gamma(ag, bg) prior over param anb and r
  bg1 <- 1
  ag1 <- 1
  bg2 <- 1
  ag2 <- 1
  
  #Beta(ab,bb)prior over q and p
  ab1 <- 2
  bb1 <- 2
  ab2 <- 2
  bb2 <- 2
  hpriorpar <- c(ag1, bg1, ag2, bg2, ab1, bb1, ab2, bb2)
  
  # Beta prior for distortion parameters (betas)
  hpriords <- c(cb, db)
  
  w <- 1  # Size of the steps for creating interval for Slice sampler (default 1)
  m <- 10 # Limit on steps for Slice sampler(default infinite)
  lo <- c(0, 0, 0, 0) # Lower bound on support of the distribution of prior parameters
  up <- c(Inf, 1, Inf, 1) # Upper bound on support of the distribution of prior parameters
  x1 <- param
  
  init <- c(anb, qnb, rnb, pnb, betas)
} else if(lam_prior == "NBD") {
  # NBD microclust prior
  samind <- c(0, 0, 0, 1, 1, 1) ## indicator vector: which parameters to slice sample
  
  anb <- 1
  qnb <- 1 - 2/N
  alpha <- 100
  rnb <- 1
  pnb <- 0.5
  param <- c(anb, qnb, alpha, rnb, pnb)
  
  #Gamma(ag1, bg1) prior over params
  ag1 <- 1
  bg1 <- 1/100
  ag2 <- 1
  bg2 <- 1
  ab1 <- 2
  bb1 <- 2
  hpriorpar <- c(ag1, bg1, ag2, bg2, ab1, bb1)
  
  # Beta prior for distortion parameters (betas)
  hpriords <- c(cb, db)
  
  w <- 1 #3,10 # Size of the steps for creating interval for Slice sampler (default 1)
  m <- 10 #30,100 # Limit on steps for Slice sampler(default infinite)
  lo <- c(0, 0, 0, 0, 0) # Lower bound on support of the distribution of prior parameters
  up <- c(Inf, 1, Inf, Inf, 1) # Upper bound on support of the distribution of prior parameters
  x1 <- param
  
  init <- c(anb, qnb, alpha, rnb, pnb, betas)
} 

# run RL ----
chaperones("RLdata500", N, lam_prior, x1, init, x, z, id, cat_fields, string_fields, 
           proportions, str_proportions, nsamples, spacing, thin, thin1, burn, 
           w, m, lo, up, hpriorpar, betas, lods, upds, hpriords, samind, rep, "results/micro_test")

