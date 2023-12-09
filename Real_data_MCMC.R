#Packages
library(RcppHMM)
library(mvtnorm)
library(bvartools)
library(mvnfast)
library(LaplacesDemon)
library(combinat)
library(DirichletReg)
library(label.switching)

#Functions
source("Marginal_likelihood_MCMC.R")
source("Data_augmented_MCMC.R")

#Faithful dataset
data("faithful")
y <- t(faithful)

#Number of iterations MCMC
#For faster computation, the number of MCMC iterations  can be reduced
iterations_m <- 300000
iterations_a <- 30000

#Number of States
states <- 2

#Scalar proposals
scalar_mu <- 0.08
scalar_sigma <- 0.14
scalar_omega <- 0.35

#Seed
set.seed(123456789)

#Marginal
M2 <- MarginalMCMC(y, iterations_m, states, scalar_mu, scalar_sigma, scalar_omega)

#Post-process with package label.switching
mcmc.pars <- array(c(t(M2$muest[1, 1:states,]), t(M2$muest[2, 1:states,]), t(M2$sigmaest[1, 1, 1:states,]), t(M2$sigmaest[2, 1, 1:states,]), t(M2$sigmaest[2, 2, 1:states,]), M2$pesest[1, 1,], M2$pesest[2, 2,]), dim = c(iterations_m, states, 6))      
run <- pra(mcmc.pars, mcmc.pars[M2$posi,,])
for(r in 1:iterations_m){
  M2$muest[,,r] <- M2$muest[, run$permutations[r,], r]
  M2$sigmaest[,,,r] <- M2$sigmaest[,, run$permutations[r,], r]
  M2$pesest[,,r] <- M2$pesest[run$permutations[r,], run$permutations[r,], r]
  M2$omegaest[,,r] <- M2$omegaest[run$permutations[r,], run$permutations[r,], r]
}

tMatrix_iniz <- M2$omegainiz/rowSums(M2$omegainiz)
mu_iniz <- M2$muiniz
sigma_iniz <- M2$variniz[,,1]

#Seed
set.seed(123456789)

#Augmented 
A2 <- AugmentedMCMC(y, iterations_a, states, "no random", 2, tMatrix_iniz, mu_iniz, sigma_iniz)

