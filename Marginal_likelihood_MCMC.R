MarginalMCMC <- function(data, iterations, states, scalar_mu, scalar_sigma, scalar_omega){
  
  y <- data
  
  #Y_t dimension
  p <- nrow(y)
  
  #Y_t sample size
  n <- ncol(y)
  
  #Log-likelihood
  if(states == 1){
    loglkd <- function(y, mu, sigma){
      log <- sum(loglik_normal((y-mu), sigma))
      return(log)
    }
  }else{
    loglkd <- function(y, mu, sigma, omega, regimes){
      Y <- y
      N <-  seq(1, regimes, 1)
      M <- nrow(mu)
      Sigma <- sigma
      tprob <- matrix(NA, regimes, regimes)
      for(j in 1:regimes){
        for(i in 1:regimes){
          tprob[j,i] <- omega[j,i]/sum(omega[j,])
        }
      }
      A <- tprob
      Mu <- mu
      Pi <- rep(1/length(N),length(N))
      HMM.cont.multi <- verifyModel(list( "Model" = "GHMM",
                                          "StateNames" = N,
                                          "A" = A,
                                          "Mu" = Mu,
                                          "Sigma" = Sigma,
                                          "Pi" = Pi))
      evaluation(HMM.cont.multi , Y , method = "f" )
    }
  }
  
  #Prior mu
  logpri_mu <- function(pr, mup, sigp){
    res <- 0
    for(i in 1:states){
      res <- res + dmvn(pr[,i], mup, sigp, log = TRUE)
    }
    return(res)
  }
  
  #Prior sigma
  logpri_sigma <- function(pr, a, lam){
    res <- 0
    for(i in 1:states){
      res <- res + dinvwishart(pr[,,i], nu = a, S = lam, log = TRUE)
    }
    return(res)
  }
  
  #Prior omega
  logpri_omega <- function(pr, a, b, c){
    res <- 0
    for(i in 1:states){
      for(j in 1:states){
        if(i==j){
          res <- res + dgamma(pr[i,j], shape = a, scale = b, log = TRUE)
        }else{
          res <- res + dgamma(pr[i,j], shape = c, scale = b, log = TRUE)
        }
      }
    }
    return(res)
  }
  
  #Jacobian sigma
  logJ_sigma <- function(sigma){
    res1 <- 1
    res2 <- 1
    for(i in 1:states){
      for(h in 1:p){
        res1 <- res1*sigma[h,h,i]
      }
    }
    for(j in 1:states){
      for(k in 1:(p-1)){
        for(l in (k+1):p){
          res2 <- res2*(sigma[k,k,j]*sigma[l,l,j]-sigma[k,l,j]^2)/(2*sqrt(sigma[k,k,j]*sigma[l,l,j]))
        }
      }
    }
    return(log(res1*res2))          
  }
  
  #Jacobian omega
  logJ_omega <- function(omega){
    res <- 1
    for(i in 1:states){
      for(j in 1:states){
        res <- res*omega[i,j]
      }
    }
    return(log(res))          
  }
  
  #Chains
  muest <- array(NA, dim = c(p, states, iterations))
  sigmaest <- array(NA, dim = c(p, p, states, iterations))
  omegaest <- array(NA, dim = c(states, states, iterations))
  
  #Prior parameters mu
  mup <- rowMeans(y)
  sigp <- matrix(0, p, p)
  #Under positive correlation assumption, if negative: sigp[si,ci] <- -(0.5*(diff(range(y[si,]))*diff(range(y[ci,])))^0.5) sigp[ci,si] <- -(0.5*(diff(range(y[si,]))*diff(range(y[ci,])))^0.5)
  for(si in 1:p){
    sigp[si,si] <- diff(range(y[si,]))
    if(si<p){
      for(ci in (si+1):p){
        sigp[si,ci] <- (0.5*(diff(range(y[si,]))*diff(range(y[ci,])))^0.5)
        sigp[ci,si] <- (0.5*(diff(range(y[si,]))*diff(range(y[ci,])))^0.5)
      }
    }
  }
  
  #Prior parameters sigma under positive correlation assumption, if negative S[al,cal] <- -half_alpha S[cal,al] <- -half_alpha
  half_alpha <- ceiling((p+1)/2) + 1
  alpha <- half_alpha*2
  S <- matrix(0, p, p)
  for(al in 1:p){
    S[al,al] <- alpha
    if(al<p){
      for(cal in (al+1):p){
        S[al,cal] <- half_alpha
        S[cal,al] <- half_alpha
      }
    }
  }
  
  #Prior parameters omega
  shape <- 1
  scale <- 1
  shape2 <- 1
  
  #Pivotal point vector
  argmax <- rep(NA, iterations)
  
  #Acceptance index
  acc <- 0
  accsig <- 0
  accom <- 0
  
  #Starting values
  muiniz <- rmvnorm(states, mean = mup, sigma = sigp)
  muiniz <- t(muiniz)
  variniz <- array(c(rep(rinvwishart(alpha,S))), dim = c(p,p,states))
  omegainiz <- matrix(c(rgamma((states*states),shape = shape, scale = scale)), states,states)
  
  #Init theta
  omu_ih <- muiniz
  osigma_ihh <- variniz
  oomega_ij <- omegainiz
  muest[,,1] <- omu_ih
  sigmaest[,,,1] <- osigma_ihh
  omegaest[,,1] <- oomega_ij
  
  #Pivotal
  if(states == 1){
    argmax[1] <- loglkd(y, muest[,,1], sigmaest[,,,1])+logpri_mu(as.matrix(muest[,,1],p,1), mup, sigp)+logpri_sigma(array(sigmaest[,,,1], dim = c(p,p,1)), alpha, S)
  }else{
    argmax[1] <- loglkd(y, muest[,,1], sigmaest[,,,1], omegaest[,,1], states)+logpri_mu(muest[,,1], mup, sigp)+logpri_sigma(sigmaest[,,,1], alpha, S)
  }
  
  #Algorithm
  time <- system.time(for(o in 2:iterations){
    nmu_ih <- matrix(NA, p, states)
    nsigma_ihh <- array(NA, dim = c(p, p, states))
    nomega_ij <- matrix(NA, states, states)
    rej <- 0
    for(i in 1:states){
      for(h in 1:p){  
        nmu_ih[h,i] <- omu_ih[h,i] + rnorm(1, 0, scalar_mu)
        nsigma_ihh[h,h,i] <- exp(log(osigma_ihh[h,h,i]) + rnorm(1, 0, scalar_sigma))
      }
    }
    for(i in 1:states){
      for(k in 1:(p-1)){
        for(l in ((k+1):p)){
          a <- sqrt(osigma_ihh[k,k,i]*osigma_ihh[l,l,i])
          an <- sqrt(nsigma_ihh[k,k,i]*nsigma_ihh[l,l,i])
          f <- log((osigma_ihh[k,l,i] + a)/(a-osigma_ihh[k,l,i])) + rnorm(1, 0, scalar_sigma)
          nsigma_ihh[k,l,i] <- (exp(f)*an -an)/(1+exp(f))
        }
      }
      nsigma_ihh[,,i][lower.tri(nsigma_ihh[,,i], diag = FALSE)] <- t(nsigma_ihh[,,i])[lower.tri(nsigma_ihh[,,i], diag = FALSE)]
    }
    for(i in 1:states){
      for(j in 1:states){
        nomega_ij[j,i] <- exp(log(oomega_ij[j,i]) + rnorm(1, 0, scalar_omega))
      }
    }
    for(i in 1:dim(nsigma_ihh)[3]){
      check <- is.positive.semidefinite(nsigma_ihh[,,i])
      if(check==TRUE){
        t <- 0 
      }else{
        t <- 1
      }
      rej <- rej + t
    }
    if(rej > 0){
      osigma_ihh <- osigma_ihh
      ratio_mu <- logpri_mu(nmu_ih,  mup, sigp) - logpri_mu(omu_ih,  mup, sigp)
      if(states == 1){
        ratio_lkd <- loglkd(y, c(nmu_ih), osigma_ihh[,,1]) - loglkd(y, c(omu_ih), osigma_ihh[,,1])
      }else{
        ratio_lkd <- loglkd(y, nmu_ih, osigma_ihh, oomega_ij, states) - loglkd(y, omu_ih, osigma_ihh, oomega_ij, states)
      }
      Amu <- ratio_lkd + ratio_mu
      
      #Update (i) if sigma is rejected
      
      #Mu
      u <- runif(1,0,1)
      if(u <= min(1,exp(Amu))){
        omu_ih <-  nmu_ih
        acc <-  acc + 1
      }else{
        omu_ih <-  omu_ih
      }
      muest[,,o] <- omu_ih
      
      #Omega
      ratio_omega <- logpri_omega(nomega_ij, shape, scale, shape2) - logpri_omega(oomega_ij, shape, scale, shape2)
      if(states == 1){
        ratio_lkd <- loglkd(y, muest[,,o], osigma_ihh[,,1]) - loglkd(y, muest[,,o], osigma_ihh[,,1])
      }else{
        ratio_lkd <- loglkd(y, muest[,,o], osigma_ihh, nomega_ij, states) - loglkd(y, muest[,,o], osigma_ihh, oomega_ij, states)
      }
      Aomega <- ratio_lkd + ratio_omega + (logJ_omega(nomega_ij) - logJ_omega(oomega_ij))
      u <- runif(1,0,1)
      if(u <= min(1,exp(Aomega))){
        oomega_ij <-  nomega_ij
        accom <-  accom + 1
      }else{
        oomega_ij <-  oomega_ij
      }
      omegaest[,,o] <- oomega_ij
    }else{
      ratio_mu <- logpri_mu(nmu_ih, mup, sigp) - logpri_mu(omu_ih, mup, sigp)
      if(states == 1){
        ratio_lkd <- loglkd(y, c(nmu_ih), osigma_ihh[,,1]) - loglkd(y, c(omu_ih), osigma_ihh[,,1])
      }else{
        ratio_lkd <- loglkd(y, nmu_ih, osigma_ihh, oomega_ij, states) - loglkd(y, omu_ih, osigma_ihh, oomega_ij, states)
      }
      Amu <- ratio_lkd + ratio_mu
      
      #Update (i) if sigma is not rejected 
      
      #Mu
      u <- runif(1,0,1)
      if(u <= min(1,exp(Amu))){
        omu_ih <-  nmu_ih
        acc <-  acc + 1
      }else{
        omu_ih <-  omu_ih
      }
      muest[,,o] <- omu_ih
      
      #Omega
      ratio_omega <- logpri_omega(nomega_ij, shape, scale, shape2) - logpri_omega(oomega_ij,shape, scale, shape2)
      if(states == 1){
        ratio_lkd <- loglkd(y, muest[,,o], osigma_ihh[,,1]) - loglkd(y, muest[,,o], osigma_ihh[,,1])
      }else{
        ratio_lkd <- loglkd(y, muest[,,o], osigma_ihh, nomega_ij, states) - loglkd(y, muest[,,o], osigma_ihh, oomega_ij, states)
      }
      Aomega <- ratio_lkd + ratio_omega + (logJ_omega(nomega_ij) - logJ_omega(oomega_ij))
      u <- runif(1,0,1)
      if(u <= min(1,exp(Aomega))){
        oomega_ij <-  nomega_ij
        accom <-  accom + 1
      }else{
        oomega_ij <-  oomega_ij
      }
      omegaest[,,o] <- oomega_ij
      
      #Sigma
      ratio_sigma <- logpri_sigma(nsigma_ihh,  alpha, S) - logpri_sigma(osigma_ihh,  alpha, S)
      if(states == 1){
        ratio_lkd <- loglkd(y, muest[,,o], nsigma_ihh[,,1]) - loglkd(y, muest[,,o], osigma_ihh[,,1])
      }else{
        ratio_lkd <- loglkd(y, muest[,,o], nsigma_ihh, omegaest[,,o], states) - loglkd(y, muest[,,o], osigma_ihh, omegaest[,,o], states)
      }
      Asigma <- ratio_lkd + ratio_sigma + (logJ_sigma(nsigma_ihh) - logJ_sigma(osigma_ihh))
      u <- runif(1,0,1)
      if(u <= min(1,exp(Asigma))){
        osigma_ihh <-  nsigma_ihh
        accsig <-  accsig + 1
      }else{
        osigma_ihh <-  osigma_ihh
      }
      sigmaest[,,,o] <- osigma_ihh
    }
    
    #Permutation
    if(states >= 2){
      ind <- sample(1:states,states)
      muest[,,o] <- omu_ih[,ind]
      omegaest[,,o] <- oomega_ij[ind,ind]
      sigmaest[,,,o] <- osigma_ihh[,,ind]
    }
    if(states == 1){
      argmax[o] <- loglkd(y, muest[,,o],sigmaest[,,,o])+logpri_mu(as.matrix(muest[,,o],p,1),mup,sigp)+logpri_sigma(array(sigmaest[,,,o],dim = c(p,p,1)),alpha,S)
    }else{
      argmax[o] <- loglkd(y, muest[,,o],sigmaest[,,,o],omegaest[,,o],states)+logpri_mu(muest[,,o],mup,sigp)+logpri_sigma(sigmaest[,,,o],alpha,S)
    }
  })
  
  #Position pivotal
  posi <- which.max(argmax)
  
  #Pi estimate
  if(states == 1){
    psest <- rep(1,iterations)
  }else{
    psest <- array(0, dim = c(states,states,iterations))
    for(i in 1:iterations){
      psest[,,i] <- omegaest[,,i]/rowSums(omegaest[,,i])
    }
  }
  res <- list(muest = muest, sigmaest = sigmaest, omegaest = omegaest, pesest = psest, rate_mu = acc/iterations, rate_sigma = accsig/iterations, rate_omega = accom/iterations, muiniz = muiniz, variniz = variniz, omegainiz = omegainiz, time = time, posi = posi, mup = mup, sigmap =sigp, alpha = alpha, S = S)
}
