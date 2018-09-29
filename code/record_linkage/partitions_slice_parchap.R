# Calculates missing, true and false links for each iteration
calcErrorOnePar <- function(z, true.links.pair, N){
    labels <- unique(z)
    if(length(labels)==N){
        est.links.pair <- list()
    }else{
        est.links <- lapply(labels, function(x) which(z==x))
        links0 <-  est.links[lapply(est.links,length) > 1]
        pairs <- lapply(links0, function(x) combn(x,2))
        upairs <- unlist(pairs)
        if(is.null(upairs)==0){
            pairs0 <- matrix(upairs,length(upairs)/2, 2, byrow=T)
            est.links.pair <- as.list(data.frame(t(pairs0)))
        }else{
            est.links.pair <- NULL
        }
    }
    # If there are no links
    if(length(est.links.pair)==0 && length(true.links.pair)==0){
        missing.links <- 0
        false.links <- 0
        true.links <- 0
    }else{
        # Correct, incorrect, and missing links
        comparison <- links.compare(est.links.pair,true.links.pair,counts.only=TRUE)
        missing.links <- comparison$missing
        true.links<-comparison$correct
        false.links <- comparison$incorrect
    }

    return(c(missing.links, true.links, false.links))
}

## Wrapper function ##
# Inputs: N - the sample size based on the dataset           #
#         Prior - Specify Prior Type: PERPS, NBNB, DP,PY,MFM #
#         param - Parameter for prior type specified         #
# Output: Parameters for probability of existing cluster (A) #
#         and new cluster (B) for reseating algorithms       #
SetParam <- function(Prior, param, N, mus = NULL, M = NULL){
    if(Prior == "PERPS") {
      # Do this if the slice sampler is on alpha and lambda #
      alpha <- param[1]
      lambda <- param[2]
      gamma <- alpha * exp(-lambda)
      A <- rep(1, N)
      B <- rep(gamma, N)
    } else if(Prior == "DP") {
    	theta <- param[1]
      A <- seq(N)
		  B <- rep(theta,N)
    } else if (Prior == "PY") {
    	theta <- param[1]
    	delta <- param[2]
    	A <- seq(N) - delta
    	B <- rep(theta, N) + seq(N) * delta
    } else if (Prior == "NBNBF"){ #new NBNB formulation as presented in nips paper
        anb <- param[1]
        qnb <- param[2]
        rnb <- param[3]
        pnb <- param[4]
        beta <- (qnb * ((1-pnb)^rnb))/(1-((1-pnb)^rnb))
        A <- (seq(N) + rnb)
        B <- (seq(N) + anb) * beta * rnb
    } else if (Prior == "NBD"){ # old NBD that didn't work
      a <- param[1]
      q <- param[2]
      A <- (seq(N) + 1)*(c((mus[-1])/mus[-(M+1)], rep(0, N-M)))
      B <- (seq(N) + a)* q * mus[1]
    }
    return(list(A=A, B=B))
}

## Wrapper function ##
# Inputs: Prior - Specify Prior Type: PERPS, NBNB, DP,PY     #
# Output: Parameter names used to save posterior samples     #
NameParamDS <- function(Prior, nfields){
    if (Prior=="DP"){
        cnames <- c("theta")#rep("delta", nfields)
        for(i in 1:nfields){
            cnames <- c(cnames, sprintf("beta_%d",i))
        }
    }
    if (Prior=="PY"){
        cnames <- c("theta", "delta")
        for(i in 1:nfields){
            cnames <- c(cnames, sprintf("beta_%d",i))
        }
    }
    if (Prior=="PERPS"){
        #cnames <- "gamma"
        cnames <- c("alpha", "lambda", "beta")
    }
    if (Prior=="NBNB"){
        #cnames <- c("a", "r", "beta")
        cnames <- c("a", "q", "r", "p")
        for(i in 1:nfields){
            cnames <- c(cnames, sprintf("beta_%d",i))
        }
    }
    if (Prior=="NBNBF"){
        #cnames <- c("a", "r", "beta")
        cnames <- c("a", "q", "r", "p")
        for(i in 1:nfields){
            cnames <- c(cnames, sprintf("beta_%d",i))
        }
    }
    if (Prior=="NBNBO"){
        #cnames <- c("a", "r", "beta")
        cnames <- c("a", "q", "r", "p")
        for(i in 1:nfields){
            cnames <- c(cnames, sprintf("beta_%d",i))
        }
    }
    if (Prior=="DMNB"){
        cnames <- c("a", "q", "alpha", "r", "p", "beta")
    }
    if (Prior=="NBD"){
        cnames <- c("a", "q", "alpha", "r", "p")
        for(i in 1:nfields){
            cnames <- c(cnames, sprintf("beta_%d",i))
        }
    }
    return(cnames)
}

## run chaperones ##
chaperones <- function(Data, N, Prior, x1, init, x, z, id, cat_fields, string_fields, 
                       proportions, str_proportions, nsamples, 
                       spacing, thin, thin1, burn, w, m, lo, up, hpriorpar,
                       betas, lods, upds, hpriords, samind, rep, out.dir){
  nfields <- ncol(x)
  thinsam <- (nsamples + burn)/thin
  thinsam1 <- nsamples/thin1

  cat(NameParamDS(Prior, nfields),"\n", file=sprintf("%s/PosParam_%s_%s_%s.txt", out.dir, Prior, Data, rep), sep=" ", append=TRUE)
  cat(c("K", "maxNk", "meanNk", "singles", "q90", "q95"), "\n", file=sprintf("%s/StatsSlice_%s_%s_%s.txt", out.dir, Prior, Data, rep), sep=" ", append=TRUE)

  lx <- loglikxSP_str(betas, as.matrix(x[, cat_fields]), as.matrix(x[, str_fields]), str_fields, cat_fields, z, proportions, str_proportions)
  kk <- 1
  ll <- 1

  meanRates <- 0
  meanK <- 0
  sumK2 <- 0

  Khat <- length(unique(z))
  tz <- table(z)
  Nk <- as.vector(tz)

  if(Prior == "NBD") {
    if(max(tz) >= M){
      M <- max(tz)
    }
    Lm <- table(factor(tz, levels=1:M))
    alpha0 = x1[3];
    r0 = x1[4];
    p0 = x1[5];
    lmu0 = lgamma(seq(M) + r0) + (seq(M))*log(p0) + r0*log(1-p0) - log(1-(1-p0)^r0) - lgamma(r0) - lgamma(seq(M) + 1)
    mu0 = exp(lmu0)
  }

  starttime <- proc.time()
  for (i in 1:(burn + nsamples)) {
    # univariate slice sampler for all parameters
    x1 <- ifelse(Prior == "NBD", 
                 unislicemNBD(x1, lx, Lm, mu0, hpriorpar, w, m, lo, up, samind), 
                 unislicem(x1, N, Khat, lx, Nk, hpriorpar, w, m, lo, up, Prior, samind))
    if(samind[length(x1) + 1] == 1) {
      betas <- unislicespb1(betas, x, z, proportions, hpriords, w, m, lods, upds, x1, N, Khat, Nk, hpriorpar, Prior)
    }

    if(Prior == "NBD") {
      #parameters of reseating algortihm #
      alpha0 = x1[3];
      r0 = x1[4];
      p0 = x1[5];
      lmu0 = lgamma(seq(M) + r0) + (seq(M))*log(p0) + r0*log(1-p0) - log(1-(1-p0)^r0)-
      lgamma(r0) - lgamma(seq(M) + 1)
      mu0 = exp(lmu0)
      alphasDM <- c(alpha0*mu0 + Lm, abs(alpha0*(1-sum(mu0))))
      mus <- rdirichlet(1, alphasDM)
      AB <- SetParam(Prior, x1, N, mus, M)
    } else {
      # parameters of reseating algortihm #
      AB <- SetParam(Prior, x1, N)
    }
    A <- AB$A
    B <- AB$B
    nclus <- length(unique(z))


    zno <- as.numeric(factor(z, labels = seq(1,length(unique(z)))))
    tclus <- nclus - max(zno)
    z1 <- Web_SamplerSP(as.matrix(x[, cat_fields]), as.matrix(x[, str_fields]), cat_fields, str_fields, 
                        zno, A, B, betas, proportions, str_proportions, 1, spacing)
    z <- as.numeric(factor(z1[1, ], labels = seq(1, length(unique(z1[1, ])))))

    if(Prior == "DP") {
      Khat <- length(unique(z))
      tz <- table(z)
      Nk <- as.vector(tz)
      if(max(tz) >= M){
          M <- max(tz)
      }
      Lm <- table(factor(tz, levels=1:M))
    }

    lx <- loglikxSP_str(betas, as.matrix(x[, cat_fields]), as.matrix(x[, str_fields]), str_fields, cat_fields, z, proportions, str_proportions)
    if(i %% thin == 0) {
      cat(c(x1, betas), "\n", file=sprintf("%s/PosParam_%s_%s_%s.txt", out.dir, Prior, Data, rep), sep=" ", append=TRUE)
      cat(lx, file=ssprintf("%s/loglik_%s_%s_%s.txt", out.dir, Prior, Data, rep),sep="\n", append=TRUE)

      missing.links <- 0
      true.links <- 0
      false.links <- 0
      for(j in 1:length(txBL)){
        if(is.null(upairs[[j]])==0){
          pairs0 <- matrix(upairs[[j]],length(upairs[[j]])/2, 2, byrow=T)
          true.links.pair <- as.list(data.frame(t(pairs0)))
          links.bl <- calcErrorOnePar(zz[[j]], true.links.pair, N)
          missing.links <- missing.links + links.bl[1]
          true.links <- true.links + links.bl[2]
          false.links <- false.links + links.bl[3]
        }
      }
      
      fnr <- missing.links/(true.links + missing.links)
      fdr <- ifelse(true.links + false.links == 0, 0, false.links/(true.links + false.links))
      rates <- c(fnr, fdr)
      cat(rates, "\n", file=sprintf("%s/rates_%s_%s_%s.txt", out.dir, Prior, Data, rep), sep=" ", append=TRUE)

      K <- length(unique(z))
      Nks <- table(z)
      maxNk <- max(Nks) # Maximum cluster size #
      meanNk <- mean(Nks) # Mean cluster size #
      singles <- sum(Nks == 1) # Number of singleton clusters #

      # Quantiles 90% and 95%
      q90 <- quantile(Nks, prob=0.9)
      q95 <- quantile(Nks, prob=0.95)
      cat(c(K, maxNk, meanNk, singles, q90, q95), "\n", file=sprintf("%s/StatsSlice_%s_%s_%s.txt", out.dir, Prior, Data, rep), sep=" ", append=TRUE)

      if(i > burn) {
          meanRates <- meanRates + rates/thinsam1
          meanK <- meanK + K/thinsam1
          sumK2 <- sumK2 + K*K
      }
    }
  }
  runtime = (proc.time() - starttime)

  #save last iteration
  write.matrix(z, file=sprintf("%s/Zetas_%s_%s_%s.txt", out.dir, Prior, Data, rep))

  post_mean_dpm <- meanK
  post_sd_dpm <- sqrt((1/(thinsam1-1))*(sumK2-thinsam1*(meanK*meanK)))
  param.to.save <- paste(init, collapse="/")
  write.csv(data.frame(Initial=param.to.save, Mean=post_mean_dpm, SD=post_sd_dpm, FNR=meanRates[1], FDR=meanRates[2], SamplerTime_Min=runtime[3]/60), sprintf("%s/LinkRates_%s_%s_%s.csv", out.dir, Prior, Data, rep), row.names=FALSE)
}
