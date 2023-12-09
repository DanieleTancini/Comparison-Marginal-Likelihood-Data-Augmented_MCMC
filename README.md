# Comparison-Marginal-Likelihood-Data-Augmented_MCMC

This repository contains the R functions used for the paper
"A comparison between marginal likelihood and data augmented MCMC algorithms for Gaussian hidden Markov models" by
- D. Tancini (University of Perugia, IT);  
- F. Bartolucci (University of Perugia, IT);
- S. Pandolfi (University of Perugia, IT).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Marginal_likelihood_MCMC.R ------> run a Marginal MCMC.

Input:
data -> D x T dataset, where D stands for dimension (e.g. 2 if bivariate) and T number of observation;
iterations -> MCMC number of iterations; 
states -> number of hidden states;
scalar mu -> sd proposal for mu;
scalar sigma -> sd proposal for sigma;
scalar omega -> sd proposal for omega.

Output:
muest -> mu chain generated, a D x K x R array, where D stands for dimension, K number of states, R number of iterations;
sigmaest -> sigma chain generated, a D x D x K x R array;
omegaest -> omega chain generated, a K x K x R array;
pesest -> transition chain generated, a K x K x R array;
rate_mu -> a scalar, acceptance rate;
rate_sigma -> a scalar, acceptance rate;
rate_omega -> a scalar, acceptance rate;
muiniz -> starting mu for the MCMC;
variniz -> starting variance-covariance for the MCMC;
omegainiz -> starting omega matrix for the MCMC;
time -> time to execute;
posi -> position required for the label switching algorithm;
mup -> mu prior;
sigmap -> sigma prior;
alpha -> alpha prior;
S -> S prior.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Data_augmented_MCM.R ------> run a Data augmented MCMC.

Input:
data -> D x T dataset, where D stands for dimension (e.g. 2 if bivariate) and T number of observation;
iterations -> MCMC number of iterations; 
states -> number of hidden states;
iniz -> MCMC initialization for the staring point, "random" if generated randomly, otherwise self initialize;
id_con -> position for the identifiability constraint;
P_iniz -> transition matrix initialization;
mu_iniz -> mu as a D x K matrix;
sig_iniz -> variance-covariance as a D x D x K array.

Output:
mu -> mu chain generated, a D x K x R array, where D stands for dimension, K number of states, R number of iterations;
sigma -> sigma chain generated, a D x D x K x R array;
Pi -> transition chain generated, a K x K x R array;
zs -> latent state as a R x T matrix;
time -> time to execute;
mup -> mu prior;
sigmap -> sigma prior;
alpha -> alpha prior;
S -> S prior.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Real_data_MCMC.R ------> script which contains a real data application. It provides both, Data augmented and Marginal likelihood MCMC.


