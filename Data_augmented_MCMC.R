AugmentedMCMC <- function(data, iterations, states, iniz, id_con, P_iniz, mu_iniz, sig_iniz){
  
  y <- data
  
  #Y_t sample size
  n <- ncol(y)
  
  #Y_t dimension
  d <- nrow(y)
  
  #Iterations 
  R <- iterations
  
  #Hidden states r=0
  z0 <-  rep(0,n)
  for(o in 1:n){
    u <-  sample(1:states,1)
    z0[o] <-  u
  }
  
  #Hyperparameters
  half_alpha <- ceiling((d+1)/2) + 1
  alpha <- half_alpha*2
  S <- matrix(0,d,d)
  for(al in 1:d){
    S[al,al] <- alpha
    if(al<d){
      for(cal in (al+1):d){
        S[al,cal] <- half_alpha
        S[cal,al] <- half_alpha
      }
    }
  }
  
  mup1 <- rowMeans(y)
  sigmap1 <- matrix(0,d,d)
  for(si in 1:d){
    sigmap1[si,si] <- diff(range(y[si,]))
    if(si<d){
      for(ci in (si+1):d){
        sigmap1[si,ci] <- (0.5*(diff(range(y[si,]))*diff(range(y[ci,])))^0.5)
        sigmap1[ci,si] <- (0.5*(diff(range(y[si,]))*diff(range(y[ci,])))^0.5)
      }
    }
  }
  
  #Theta_0
  iniz <- iniz
  if(iniz == "random"){
    Pi1 <- matrix(1/states, states,states)
    muiniz <- rmvnorm(states, mean = mup1, sigma = sigmap1)
    muiniz <- t(muiniz) 
    variniz <- rinvwishart(alpha,S)
  }else{
    Pi1 <- P_iniz
    muiniz <- mu_iniz
    variniz <- sig_iniz
  }
  z <- z0
  T <- length(z)
  
  #Init theta
  mu <- array(0,dim = c(d,states,R))
  sigma <- array(0,dim = c(d,d,states,R))
  Pi <- array(0, dim = c(states,states,R))
  Pi[,,1] <- Pi1
  mu[,,1] <- muiniz
  sigma[,,,1] <- variniz
  muf <- mu[,,1]
  sigmaf <- sigma[,,,1]
  Pif <- Pi1
  zs <-  matrix(0,R,T)
  zs[1,] <- z
  
  #Exact calculation full
  solvesigmap1 <- solve(sigmap1)
  solvesigmap1mup1 <- solve(sigmap1)%*%mup1
  
  #Algorithm
  dens <-  matrix(0,T,states)
  
  time <- system.time(for(r in 2:R){
    if(states == 1){
      Pi[,,r] <- 1
      zs[r,] <- 1
      mu[,1,r] <- rmvnorm(1,solve(solvesigmap1+T*solve(sigma[,,1,r-1]))%*%(solvesigmap1mup1 + solve(sigma[,,1,r-1])%*%rowSums(y)),solve((solvesigmap1+T*solve(sigma[,,1,r-1]))))
      sigma[,,1,r] <- rinvwishart(alpha + T,S + (y - mu[,1,r])%*%t(y - mu[,1,r]))
    }else{
      for(i in 1:states){
        dens[,i] <- exp(loglik_normal((y-muf[,i]), sigmaf[,,i]))
      }
      for(t in T:1){
        if(t==T){
          pz <- Pif[z[t-1],]*dens[t,]
          Pz <- pz/sum(pz)
          z[t] = which(rmultinom(1,1,Pz)==1)
        }else if(1<t & t<T){
          pz <- Pif[z[t-1],]*dens[t,]*Pif[,z[t+1]]
          Pz <- pz/sum(pz)
          z[t] = which(rmultinom(1,1,Pz)==1)
        }else{
          pz <- Pif[,z[t+1]]*dens[t,]
          Pz <- pz/sum(pz)
          z[t] = which(rmultinom(1,1,Pz)==1)
        }
      }
      aa <- as.vector(table(factor(z[1:length(z)-1], levels = 1:states), factor(z[2:length(z)], levels = 1:states)))
      for(u in 1:states){
        dir <- seq(u,(states^2),states)
        Pi[u,,r] <- rdirichlet(1, c(1+aa[dir]))
        ns <- ifelse(z==u,1,0)
        yn <- y*matrix(rep(ns,d), ncol = T, byrow = TRUE)
        mu[,u,r] <- rmvnorm(1,solve(solvesigmap1+sum(ns)*solve(sigma[,,u,r-1]))%*%(solvesigmap1mup1 + solve(sigma[,,u,r-1])%*%rowSums(yn)),solve((solvesigmap1+sum(ns)*solve(sigma[,,u,r-1]))))
        sum <- colSums(yn)
        pos <- ifelse(sum==0,1,0)
        posiz <- which(pos==1,arr.ind = T)
        ynn <- as.matrix(yn[,-c(posiz)])
        if(dim(ynn)[2]==0){
          Sn <- S
        }else{
          Sa <- array((ynn-mu[,u,r]), dim = c(d,1,ncol(ynn)))
          St <- matrix(rowSums(apply(Sa, 3, function(x) x%*%t(x))), nrow = d)
          Sn <- S + St
        }
        an <- alpha + sum(ns)
        sigma[,,u,r] <- rinvwishart(an,Sn)
      }
      
      #Identif constraint
      norm_e <- mu[id_con,,r]
      ind <-  order(norm_e)
      if(any(ind!=(1:states))){
        mu[,,r] <-  mu[,ind,r]
        sigma[,,,r] <-  sigma[,,ind,r]
        Pi[,,r] <-  Pi[ind,ind,r]
        z1 <- z
        old <- c(1:states)
        new <- ind
        z <- new[match(z1,old)]
      }
      muf <- mu[,,r]
      sigmaf <- sigma[,,,r]
      Pif <- Pi[,,r]
      zs[r,] <- z
    }
  })
  res <- list(mu = mu, sigma = sigma, Pi = Pi, zs = zs, time = time, alpha = alpha, S = S, mup = mup1, sigmap = sigmap1)
}
