rm(list = ls())
# Source all the required libraries
library(MixMatrix)
library(MCMCpack)
library(telefit)
library(NPVecchia)
library(Rcpp)
library(fields) # Contains rdist function
library(adaptMCMC)
library(nloptr)
library(tidyverse)
library(mvtnorm)
library(BDgraph)

# Quiet Function 
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

dmatrix_norm <- function(X, U, V){
  n = nrow(X); p = ncol(X)
  
  first = - 0.5 * sum(diag(solve(V, t(X)) %*% solve(U, X)))
  second = -0.5 * n * p * log(2*pi)
  third = - 0.5 * n * determinant(V, logarithm = TRUE)$modulus[1]
  forth = - 0.5 * p * determinant(U, logarithm = TRUE)$modulus[1]
  return(first + second + third + forth)
}  

# number of locations, dimensions, samples
n_locs1 <- 100
n_locs2 <- 150
n_locs3 <- 200

n_locs <- c(n_locs1, n_locs2, n_locs3)

d <- 2

num_reps <- 20

# random locations and data
locs1 <- matrix(runif(d*n_locs1, min = 0, max = 1), nc = d)
locs2 <- matrix(runif(d*n_locs2, min = 0, max = 1), nc = d)
locs3 <- matrix(runif(d*n_locs3, min = 0, max = 1), nc = d)

# Column covariances
Sigma1.true = maternCov(d = fields::rdist(locs1), range = 2, smoothness = 0.25, scale = 1)
Sigma2.true = maternCov(d = fields::rdist(locs2), range = 2, smoothness = 0.5, scale = 1.5)
Sigma3.true = maternCov(d = fields::rdist(locs3), range = 2, smoothness = 0.5, scale = 2)

# Row covariances
##############################
# AR correlation structure
##############################
r = 0.5
Psi = matrix(0, num_reps, num_reps) # Define the true scale matrix
# The true scale matrix is generated as an AR type matrix
for(i in 1:num_reps){
  for(j in 1:num_reps){
    Psi[i, j] = r^(abs(i - j))
  }
}

Lambda.true = as.matrix(forceSymmetric(riwish(v = num_reps, S = Psi)))

# Generate the data from a matrix normal
dataa1 <- MixMatrix::rmatrixnorm(n = 1, mean = matrix(0, nrow = num_reps, ncol = n_locs1), U = Lambda.true, V = Sigma1.true)

dataa2 <- MixMatrix::rmatrixnorm(n = 1, mean = matrix(0, nrow = num_reps, ncol = n_locs2), U = Lambda.true, V = Sigma2.true)

dataa3 <- MixMatrix::rmatrixnorm(n = 1, mean = matrix(0, nrow = num_reps, ncol = n_locs3), U = Lambda.true, V = Sigma3.true)

ll.true = dmatrix_norm(X = dataa1, U = Lambda.true, V = Sigma1.true) + dmatrix_norm(X = dataa2, U = Lambda.true, V = Sigma2.true) + dmatrix_norm(X = dataa3, U = Lambda.true, V = Sigma3.true)

##### PRE-PROCESSING THE DATA
# Compute maximin ordering
order1 <- order_maximin_dist(fields::rdist(locs1))
order2 <- order_maximin_dist(fields::rdist(locs2))
order3 <- order_maximin_dist(fields::rdist(locs3))

# Reorder data and location by maximin ordering
dataa1 <- dataa1[, order1]; dataa2 <- dataa2[, order2]; dataa3 <- dataa3[, order3]
locs1 <- locs1[order1, ]; locs2 <- locs2[order2, ]; locs3 <- locs3[order3, ]

# Find the Euclidean neighbors
nearest_neighbors1 <- find_nn_dist(fields::rdist(locs1), n_locs1)
nearest_neighbors2 <- find_nn_dist(fields::rdist(locs2), n_locs2)
nearest_neighbors3 <- find_nn_dist(fields::rdist(locs3), n_locs3)

n_iterations = 1000 # Number of iterations


source("functions_replicates.R") 
#####################################################################################################
# My method
#####################################################################################################
dataa <- list(dataa1, dataa2, dataa3)
nearest_neighbors <- list(nearest_neighbors1, nearest_neighbors2, nearest_neighbors3)
m = 10 # Number of nearest neighbors 
N = nrow(dataa1) # Define the number of replicates

nu = N # Degrees of freedom for the IW prior

Psi = diag(1, num_reps) # Set scaling parameter for IW prior
Lambda = as.matrix(forceSymmetric(riwish(v = num_reps, S = Psi))) # Starting value for the IW sampling (taken to be equal to true value)

# Define lists to store the samples
U_samples <- replicate(length(dataa), list(), simplify = FALSE)
d_samples <- replicate(length(dataa), list(), simplify = FALSE)

post <- list()

theta_samples <- matrix(NA, nrow = n_iterations, ncol = 3)

Lambda_samples <- list()

#A scale matrix similar to this worked well for our application
scale_mat <- matrix(c(0.05, -0.04, 0, -0.04, 0.05, 0, 0, 0, 0.01), nc = 3)

#number of samples
nruns <- 2
log_likelihood <- list()
# Lambda = diag(1, nrow = num_reps)
init_theta = c(1,-1,0)

start_time <- Sys.time()
for(i in 1:n_iterations){
  # Printing the iterations
  if(i == 1){
    cat(paste0("Iteration: ", i, "\n"))
  }
  if(i %% floor((5/100)*(n_iterations + 1)) == 0) {
    cat(paste0("Iteration: ", i, "\n"))
  }
  # Find the optimal thetas as a function of Lambda from the previous iteration
  thetas_new <- quiet(adaptMCMC::MCMC(minus_loglikeli_my_new, datumT = dataa, 
                                      NNarrayT = nearest_neighbors,
                                      Lambda = Lambda,
                                      m = m,
                                      init = init_theta, negativ = FALSE,
                                      scale = scale_mat, adapt = TRUE, 
                                      acc.rate = 0.234, n = nruns, 
                                      showProgressBar = FALSE)$samples[nruns, ])
  
  theta_samples[i, ] = thetas_new
  for(r in 1:length(dataa)){
    # Convert the hyper-parameters thetas to priors which is a function of Lambda from previous iteration
    priors = thetas_to_priors_my(thetas_new, n = nrow(nearest_neighbors[[r]]), m = m)
    # get initial posterior sample
    post[[r]] = get_posts_my(datum = dataa[[r]], priors = priors, Lambda = Lambda, NNarray = nearest_neighbors[[r]])
    
    # Sample posterior sparse matrix U as a function of Lambda from previous iteration
    U_samples[[r]][[i]] = samp_posts_my_new(post[[r]], NNarray = nearest_neighbors[[r]], bayesian = TRUE)
    # Store the diagonal elements of U
    d_samples[[r]][[i]] = diag(U_samples[[r]][[i]])
  }
  
  
  SS = matrix(0, nrow = N, ncol = N) # Define the sample sum of square matrix
  for(r in 1:length(dataa)){
    for(ind in 1:n_locs[r]){
      if(ind == 1){
        mean.normal = rep(0, num_reps)
        d = 1/((d_samples[[r]][[i]][ind])^2)
        SS = SS + tcrossprod(dataa[[r]][, ind])/d
      }else{
        X = as.matrix(post[[r]][[5]][[ind]])
        d = 1/((d_samples[[r]][[i]][ind])^2)
        u =  sqrt(d) * U_samples[[r]][[i]][as.numeric(post[[r]][[6]][[ind]]), ind]  
        mean.normal = as.numeric(X %*% u)
        SS = SS + (tcrossprod(dataa[[r]][, ind] - mean.normal))/d
      }
    }
  }
  
  # Draw sample from IW i.e. from full conditional [Lambda|-]
  Lambda = riwish(v = nu + sum(n_locs), S = as.matrix(forceSymmetric(Psi + SS)))
  # Store Lambda as sample
  Lambda_samples[[i]] = as.matrix(forceSymmetric(Lambda))
  # log_like <- list(list(), list())
  log_like <- replicate(length(dataa), list(), simplify = FALSE)
  for(r in 1:length(dataa)){
    for(ind in 1:n_locs[[r]]){
      if(ind == 1){
        mean.normal = rep(0, num_reps)
        d = 1/((d_samples[[r]][[i]][ind])^2)
      }else{
        X = as.matrix(post[[r]][[5]][[ind]])
        d = 1/((d_samples[[r]][[i]][ind])^2)
        u =  sqrt(d) * U_samples[[r]][[i]][as.numeric(post[[r]][[6]][[ind]]), ind]  
        mean.normal = as.numeric(X %*% u)
        
      }
      
      log_like[[r]][[ind]] = dmvnorm(x = dataa[[r]][, ind], mean = mean.normal, sigma = d * Lambda_samples[[i]], log = TRUE)
    }
    log_likelihood[[i]] = sum(unlist(log_like))
  }
}# End of Gibbs Sampling
end_time <- Sys.time()
end_time - start_time

U1_post = Matrix(0, nrow = n_locs1, ncol = n_locs1, sparse = TRUE)
U2_post = Matrix(0, nrow = n_locs2, ncol = n_locs2, sparse = TRUE)
U3_post = Matrix(0, nrow = n_locs3, ncol = n_locs3, sparse = TRUE)
Lambda_post = matrix(0, nrow = num_reps, ncol = num_reps)
theta_post = 0

burn = 100 # Number of burn-in

thin = 1

samples <- seq(from = (burn + 1), to = n_iterations, by = thin)

par(mfrow = c(1,1))
plot(unlist(log_likelihood[samples]), type = "l", ylab = "LL", xlab = "Iteration", main = "Plot of log-likelihood")
abline(h = ll.true, col = "red")

for(i in samples){
  U1_post = U1_post + U_samples[[1]][[i]]  # Sum of all matrices from Posterior of U1
  U2_post = U2_post + U_samples[[2]][[i]]  # Sum of all matrices from Posterior of U2
  U3_post = U3_post + U_samples[[3]][[i]]  # Sum of all matrices from Posterior of U3
  Lambda_post = Lambda_post + Lambda_samples[[i]] # Sum of all matrices from Posterior of Lambda
}

U1_post = U1_post/length(samples) # Posterior Mean as an estimate of U1
U2_post = U2_post/length(samples) # Posterior Mean as an estimate of U2
U3_post = U3_post/length(samples) # Posterior Mean as an estimate of U2
Lambda_post = Lambda_post/length(samples) # Posterior Mean as an estimate of Lambda
theta_post = apply(theta_samples, MARGIN = 2, FUN = mean)

Sigma1_post = solve(t(U1_post),  solve(U1_post)) # Posterior estimate of Sigma1
Sigma2_post = solve(t(U2_post),  solve(U2_post)) # Posterior estimate of Sigma2
Sigma3_post = solve(t(U3_post),  solve(U3_post)) # Posterior estimate of Sigma2


# Calculate \hat{\Sigma} \Sigma^{-1}
prod1 = Matrix(cov2cor(Sigma1.true) %*% solve(cov2cor(Sigma1_post)))
prod2 = Matrix(cov2cor(Sigma2.true) %*% solve(cov2cor(Sigma2_post)))
prod3 = Matrix(cov2cor(Sigma3.true) %*% solve(cov2cor(Sigma3_post)))

# Calculate the KL divergence of the column covariance matrix Sigma 
# NOTE: determinate is calculated in the log scale
KL1 = log(0.5 * (sum(diag(prod1)) - determinant(prod1, log = TRUE)$modulus[1] - ncol(Sigma1.true))) 
KL2 = log(0.5 * (sum(diag(prod2)) - determinant(prod2, log = TRUE)$modulus[1] - ncol(Sigma2.true))) 
KL3 = log(0.5 * (sum(diag(prod3)) - determinant(prod3, log = TRUE)$modulus[1] - ncol(Sigma3.true))) 

# Calculate the KL divergence for Matrix Normal 
V1 = cov2cor(Sigma1.true); V2 = cov2cor(Sigma1_post); U1 = cov2cor(Lambda.true); U2 = cov2cor(Lambda_post)  # Look at notes for definitions of U1, U2, V1, V2
# KL divergence calculation. All log(det) are calculated using determinant( , log = TRUE) and appropriate ratios are replaced by differences. Formula is complicated. Be sure to double check
KL_Matrix_Normal1 = log((sum(diag(solve(V2, V1))) * sum(diag(solve(U2, U1))) - (num_reps*(determinant(V1, log = TRUE)$modulus - determinant(V2, log = TRUE)$modulus))[1] - (n_locs1 * (determinant(U1, log = TRUE)$modulus - determinant(U2, log = TRUE)$modulus))[1] - (n_locs1 * num_reps))/2)


V1 = cov2cor(Sigma2.true); V2 = cov2cor(Sigma2_post); U1 = cov2cor(Lambda.true); U2 = cov2cor(Lambda_post)  # Look at notes for definitions of U1, U2, V1, V2
# KL divergence calculation. All log(det) are calculated using determinant( , log = TRUE) and appropriate ratios are replaced by differences. Formula is complicated. Be sure to double check
KL_Matrix_Normal2 = log((sum(diag(solve(V2, V1))) * sum(diag(solve(U2, U1))) - (num_reps*(determinant(V1, log = TRUE)$modulus - determinant(V2, log = TRUE)$modulus))[1] - (n_locs2 * (determinant(U1, log = TRUE)$modulus - determinant(U2, log = TRUE)$modulus))[1] - (n_locs2 * num_reps))/2)


V1 = cov2cor(Sigma3.true); V2 = cov2cor(Sigma3_post); U1 = cov2cor(Lambda.true); U2 = cov2cor(Lambda_post)  # Look at notes for definitions of U1, U2, V1, V2
KL_Matrix_Normal3 = log((sum(diag(solve(V2, V1))) * sum(diag(solve(U2, U1))) - (num_reps*(determinant(V1, log = TRUE)$modulus - determinant(V2, log = TRUE)$modulus))[1] - (n_locs3 * (determinant(U1, log = TRUE)$modulus - determinant(U2, log = TRUE)$modulus))[1] - (n_locs3 * num_reps))/2)


norm_Sigma1_post = norm(cov2cor(Sigma1_post) - cov2cor(Sigma1.true), type = "F")/norm(cov2cor(Sigma1.true), type = "F") 
norm_Sigma2_post = norm(cov2cor(Sigma2_post) - cov2cor(Sigma2.true), type = "F")/norm(cov2cor(Sigma2.true), type = "F") 
norm_Sigma3_post = norm(cov2cor(Sigma3_post) - cov2cor(Sigma3.true), type = "F")/norm(cov2cor(Sigma3.true), type = "F") 
norm_Lambda_post = norm(cov2cor(Lambda_post) - cov2cor(Lambda.true), type = "F")/norm(cov2cor(Lambda.true), type = "F") 


KL.estimates = list(KL1 = KL1, KL2 = KL2, KL3 = KL3,
                    KL_Matrix_Normal1 = KL_Matrix_Normal1,
                    KL_Matrix_Normal2 = KL_Matrix_Normal2,
                    KL_Matrix_Normal3 = KL_Matrix_Normal3)


norm_comparisons_no_constraint <- list(norm_Sigma1_post = norm_Sigma1_post,
                                       norm_Sigma2_post = norm_Sigma2_post,
                                       norm_Sigma3_post = norm_Sigma3_post,
                                       norm_Lambda_post = norm_Lambda_post
)



