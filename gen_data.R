# Source all the required libraries
if(!require(MixMatrix))install.packages("MixMatrix"); suppressPackageStartupMessages(library(MixMatrix))
if(!require(MCMCpack))install.packages("MCMCpack"); suppressPackageStartupMessages(library(MCMCpack))
if(!require(telefit))install.packages("telefit"); suppressPackageStartupMessages(library(telefit))
if(!require(NPVecchia))install.packages("NPVecchia"); suppressPackageStartupMessages(library(NPVecchia))
if(!require(Rcpp))install.packages("Rcpp"); suppressPackageStartupMessages(library(Rcpp))
if(!require(fields))install.packages("fields"); suppressPackageStartupMessages(library(fields)) # Contains rdist function
if(!require(adaptMCMC))install.packages("adaptMCMC"); suppressPackageStartupMessages(library(adaptMCMC))
if(!require(nloptr))install.packages("nloptr"); suppressPackageStartupMessages(library(nloptr))
if(!require(tidyverse))install.packages("tidyverse"); suppressPackageStartupMessages(library(tidyverse))
if(!require(mvtnorm))install.packages("mvtnorm"); suppressPackageStartupMessages(library(mvtnorm))
if(!require(BDgraph))install.packages("BDgraph"); suppressPackageStartupMessages(library(BDgraph))
suppressPackageStartupMessages(source("functions.R"))

# Function to simulate matrix-variate spatial data 
# @n_locs: number of spatial locations.
# @d: dimension of spatial locations. Defaults to 2. Usually 2 or 3.
# @range: range of the Matern kernel used to simulate the true column covariance matrix.
# @smoothness: smoothness parameter of the Matern kernel used to simulate the true column covariance matrix.
# @scale: marginal variance of the Matern kernel used to simulate the true column covariance matrix.
# @r: value of rho to simulate the scale matrix corresponding to the row covariance matrix.
# @type: one of "AR", "Band", "Equi". Denotes the type of row scale structure. Either AR correlation, Banded correlation or Equi correlation structure used to simulate the scale matrix corresponding to the row covariance matrix.

generate_data <- function(n_locs, d = 2, num_reps, range = 1, smoothness = 0.5, scale = 1, r = 0.5, type = "AR"){
# random locations and data
locs <- matrix(runif(d*n_locs, min = 0, max = 1), nc = d)

# Column covariances
Sigma.true = maternCov(d = fields::rdist(locs), range = range, smoothness = smoothness, scale = scale)
if(type == "AR"){
# Row covariances
##############################
# AR correlation structure
##############################
  Psi = matrix(0, num_reps, num_reps) # Define the true scale matrix
  # The true scale matrix is generated as an AR type matrix
  for(i in 1:num_reps){
    for(j in 1:num_reps){
      Psi[i, j] = r^(abs(i - j))
    }
  }
  Lambda.true = as.matrix(forceSymmetric(riwish(v = num_reps, S = Psi)))
}
##############################
# Banded correlation structure
##############################
if(type == "Band"){
  Psi = matrix(0, num_reps, num_reps) # Define the true scale matrix
  for(i in 2:(num_reps - 1)){
    for(j in (i - 1):(i + 1)){
      Psi[i, j] = r
    }
  }
  
  Psi[1, 2] = r; Psi[num_reps, (num_reps-1)] = r
  diag(Psi) = 1
  
  Lambda.true = as.matrix(forceSymmetric(riwish(v = num_reps, S = Psi)))
}
#############################
## Equicorrelation structure
############################
if(type == "Equi"){
  Psi = matrix(r, nrow = num_reps, ncol = num_reps)
  diag(Psi) = 1
  # True covariance for rows
  Lambda.true = as.matrix(forceSymmetric(riwish(v = num_reps, S = Psi)))
}

# Generate the data from a matrix normal
dataa <- MixMatrix::rmatrixnorm(n = 1, mean = matrix(0, nrow = num_reps, ncol = n_locs), U = Lambda.true, V = Sigma.true)
# Calculate the true log-likelihood values
ll.true = dmatrix_norm(X = dataa, U = Lambda.true, V = Sigma.true)
ll.true.cor = dmatrix_norm(X = dataa, U = cov2cor(Lambda.true), V = cov2cor(Sigma.true))

# Return the matrix-variate data, true row and column covariances, spatial locations and log-likelihoods based on row and column covariance and correlation matrices respectively
return(list(data = dataa, Sigma.true = Sigma.true, Lambda.true = Lambda.true, locs = locs, ll.true = ll.true, ll.true.cor = ll.true.cor))
}

