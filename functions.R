###############################################################################
# Functions in order (Auxillary functions)
###############################################################################
# Function to calculate the maximin ordering
# @d is distance matrix
order_maximin_dist_my <- function(d){
  ## number of locations to order
  n = nrow(d)
  
  ## initialize ordering
  ord = numeric(length = n)
  #get first location (basically the most central point)
  ord[1] = which.min(rowSums(d))
  
  ## compute maxmin ordering
  #Get current column
  mint <- d[,ord[1]]
  for(i in 2:n){
    #Get the maximum minimum distance
    a <- which.max(mint)
    #Add it as the next point in the ordering
    ord[i] <- a
    #Update the minimums. This works because the current points distance from itself
    #is zero, so it gets ignored when taking the maximum in the future (as all distances
    #must be >= 0 to be a valid distance).
    mint <- pmin(mint, d[,a])
  }
  #returns the ordering
  return(ord)
}

# Function to calculate the nearest neighbours
# @d: distance matrix
# @m: the number of locations
find_nn_dist_my <- function(d, m){
  n = nrow(d)
  
  ## find ordered NN
  #initialize
  NN = matrix(NA, n, m)
  for(i in 2:n){
    # if((i %% 1000)==0) print(i)
    #get number of neighbors, m if that many previous points
    k = min(i - 1, m)
    NN[i, 1:k]=order(d[i, 1:(i-1)])[1:k]
  }
  
  return(NN)
}


# For our method we need to fix m (number of neighbours) rather than estimate it from the data
# Removing the threshold parameter and calculation of m. If deemed to be not a problem
# put back threshold

# Function to convert the hyper-parameter theta to hyper-parameters of NIG prior
# @theas: vector of hyper-parameters
# @n: number of locations
# @m: number of neighbours
thetas_to_priors_my <- function(thetas, n, m) {
  # # warn if thetas are not in safe range
  # if(any(thetas > 4) || any(thetas < -6)){
  #   warning("A theta being too large/small will probably cause numerical issues.")
  # } 
  # Inverse-gamma scale prior parameter vector (prior on variances)
  b <- 5 * exp(thetas[[1]]) * (1 - exp(-exp(thetas[[2]])/sqrt(0:(n - 1))))
  # Inverse-gamma shape prior parameter vector (prior on variances)
  a <- rep(6, n)
  # temporary vector for determining number of neighbors
  tempor <- exp(-exp(thetas[[3]]) * (1:m))
  
  # Force at least 2 neighbors
  if (is.na(m) | m < 2) {
    m <- 2
  }
  # Create the prior on the coefficient variances
  g <- matrix(tempor[1:m], ncol = m, nrow = n, byrow = T)
  # Divide by mean of the IG prior for simplicity of derivations
  g <- g/(b/(a - 1))
  return(list(a, b, g))
}

# For our method we need to fix m (number of neighbours) rather than estimate it from the data
# Removing the threshold parameter and calculation of m. If deemed to be not a problem
# put back threshold

### NOTE:
# The new method forces the change to the marginal log-likelihood of theta
# Function to calculate the marginal log-likelihood
# @thetas: vector of hyper-parameters
# @datum: matrix-variate data
# @NNarray: array of nearest neighbours
# @m: number of nearest neighbours
# @Lambda: row covariance matrix
minus_loglikeli_my_new <- function(thetas, datum, NNarray, m, Lambda, negativ = TRUE) {
  # # warn if thetas are not in safe range
  # if(any(thetas > 4) || any(thetas < -6)){
  #   warning("A theta being too large/small will probably cause numerical issues.")
  # }
  # make sure datum and NNarray are both matrices
  if(!(is.matrix(datum) && is.matrix(NNarray))) {
    stop("The data and NNarray must both be matrices")
  }
  # make sure NNarray is of correct size
  if(ncol(datum) != nrow(NNarray)){
    stop(paste("The number of locations (", ncol(datum), ") must equal the number of 
               rows of the neighbor matrix but given (", nrow(NNarray), ")", sep=""))
  }
  if(ncol(NNarray) < 2) {
    stop("At least 2 neighbors are required (2 or more columns in NNarray)")
  }
  # check if NNarray is an integer matrix
  if(! is.integer(NNarray)){
    warning("NNarray should consist of only integers (and/or NAs)!")
  }
  
  # get n, N
  n <- nrow(NNarray) # This corresponds to the number of spatial locations
  N <- nrow(datum) # This corresponds to the number of replicates
  # get alpha posterior
  a_post <- 6 + N/2
  # get priors
  # Removing this and fixing m
  pr <- thetas_to_priors_my(thetas, n, m = m)
  # get m
  # make sure m is not greater than max number of neighbors set
  # (from NNarray creation)
  m <- min(ncol(pr[[3]]), ncol(NNarray))
  # get needed priors from the list of priors
  b <- pr[[2]] # Prior for beta. Here beta1 should not be used
  g <- pr[[3]] # Prior for V_i
  
  # get the first element for the log-likelihood
  loglikelihood <- 6 * log(b[1]) - a_post * log(b[1] + crossprod(datum[, 1], solve(Lambda, datum[, 1]))/2)
  # Since in the integrated likelihood the index is from 2,..,n and index i=1 does not appear. We remove the previous case of considering the log likelihood from i = 1 and then in a loop run from i = 2 to n  
  for (i in 2:n) {
    # get nearest neighbors and how many there are as nn
    gind <- na.omit(NNarray[i, 1:m])
    nn <- length(gind)
    # set up i'th regression with basic notation yi ~ xi
    xi <- -datum[, gind] 
    yi <- datum[, i]
    # Get inverse of G posterior
    Ginv <- crossprod(xi, solve(Lambda, xi)) + diag(g[i,1:nn]^(-1), nrow = nn)
    # Take the Cholesky of Ginv
    Ginv_chol <- chol(Ginv)
    # Try to get muhat directly, but if Ginv_chol is not invertible, use 
    # the generalized inverse
    muhat <- tryCatch({
      solve(Ginv_chol, solve(t(Ginv_chol), crossprod(xi, solve(Lambda, yi))))
    }, error = function(e) {
      # Replace NaN values in Ginv by 0
      Ginv[!is.finite(Ginv)] <- 0
      # pracma::pinv(Ginv) %*% crossprod(xi, solve(Lambda, yi))
      pracma::pinv(Ginv) %*% crossprod(xi, solve(Lambda, yi))  # Replace inverse of Lambda by g-inverse
    })
    # get the posterior of b (IG scale parameter)
    b_post <- b[i] + (crossprod(yi, solve(Lambda, yi)) - t(muhat) %*% Ginv %*% muhat)/2
    # b_post <- b[i] + (crossprod(yi, pracma::pinv(Lambda) %*% yi) - crossprod(muhat, Ginv %*% muhat))/2  # Replace inverse of Lambda by g-inverse
    # calculate the determinant term of the integrated likelihood
    log_det <- -0.5 * (2 * sum(log(na.omit(diag(Ginv_chol)))) + (sum(log(g[i, 1:nn]))))
    # calculate the term based on the ratio of IG parameters
    log_ig <- 6 * log(b[i]) - a_post * log(b_post)
    # Add these values to the log integrated likelihood
    loglikelihood <- loglikelihood + log_det + log_ig
  }
  # Return the negative of the log integrated likelihood
  loglikelihood <- ifelse(negativ, -1, 1)*c(loglikelihood)
  return(loglikelihood)
}

# Function to calculate the form of the posterior parameters (from the full conditional distributions)
# @datum: matrix-variate data
# @priors: list of hyper-parameters of the NIG prior (this corresponds to the output from the function thetas_to_priors_my)
# @Lambda: row covariance matrix
# @NNarray: array of nearest neighbours

get_posts_my <- function(datum, priors, Lambda, NNarray) {
  # check if priors is of correct length
  if(length(priors) != 3) {
    stop("Priors should be a list of length 3!")
  }
  # get b, g priors from the list
  b <- priors[[2]]
  g <- priors[[3]]
  # get n, N, m
  n <- ncol(datum) # Number of spatial locations
  N <- nrow(datum) # Number of replicates
  m <- min(ncol(g), ncol(NNarray)) # Force m to be the minimum of ncol(g) and ncol(NNarray)
  
  # make sure datum and NNarray are both matrices
  if(!(is.matrix(datum) && is.matrix(NNarray))) {
    stop("The data and NNarray must both be matrices")
  }
  # make sure NNarray is of correct size
  if(ncol(datum) != nrow(NNarray)){
    stop(paste("The number of locations (", ncol(datum), ") must equal the number of 
               rows of the neighbor matrix but given (", nrow(NNarray), ")", sep=""))
  }
  if(ncol(NNarray) < 2) {
    stop("At least 2 neighbors are required (2 or more columns in NNarray)")
  }
  # check if NNarray is an integer matrix
  if(! is.integer(NNarray)){
    warning("NNarray should consist of only integers (and/or NAs)!")
  }
  # check if a, b, g are correct sizes
  if(length(priors[[1]]) != n || length(b) != n || nrow(g) != n || ncol(g) < 2){
    stop("Please use priors created by thetas_to_priors for convenience. 
          The current priors are of the wrong size.")
  }
  
  # Create vectors/matrices/arrays to hold the posteriors
  a_post <- rep(0, n)
  b_post <- rep(0, n)
  muhat_post <- matrix(NA, nrow = n, ncol = m)
  G_post <- array(NA, dim = c(m, m, n))
  # Get posterior of a
  a_post <- priors[[1]] + N/2
  # Get first element of posterior of b
  b_post[1] <- b[1] + crossprod(datum[, 1], solve(Lambda, datum[, 1]))/2
  # Define a list to store the X_i 's
  X <- list()
  X[[1]] = rep(0, N)  # By definition, X_1 = 0
  indices = list() # Define a list to store the indices of non-zero off-diagonal elements of U
  indices[[1]] = 0 # By definition, the first element is 0
  for (i in 2:n) {
    # set up the regression as column i ~ nearest neighbors
    gind <- na.omit(NNarray[i, 1:m])
    indices[[i]] = gind # The index of non-zero off-diagonal elements
    nn <- length(gind) 
    xi <- -datum[, gind] 
    X[[i]] = xi  # Store X_i 's
    yi <- datum[, i]
    # Get inverse of posterior of G (coefficient variance)
    Ginv <- crossprod(xi, solve(Lambda, xi)) + diag(g[i, 1:nn]^(-1), nrow = nn)
    # take Cholesky
    Ginv_chol <- chol(Ginv)
    # G <- ginv(Ginv); muhat <- G%*%t(xi)%*%yi
    # Try to solve for muhat directly which occasionally has
    # numerical issues. If so, use generalized inverse of Ginv
    muhat <- tryCatch({
      solve(Ginv_chol, solve(t(Ginv_chol), crossprod(xi, solve(Lambda, yi))))
    }, error = function(e) {
      # Replace NaN values in Ginv by 0
      Ginv[!is.finite(Ginv)] <- 0
      pracma::pinv(Ginv) %*% crossprod(xi, solve(Lambda, yi))
    })
    # Fill in full posteriors with the elements from i'th regression
    muhat_post[i, 1:nn] <- muhat
    G_post[1:nn, 1:nn, i] <- Ginv_chol
    b_post[i] <- b[i] + (crossprod(yi, solve(Lambda, yi)) - t(muhat) %*% Ginv %*% muhat)/2
  }
  #Return list of posterior parameters
  return(list(a_post, b_post, muhat_post, G_post, X, indices))
}


# The function to sample from the full conditional distributions of U and d
# @datum: matrix-variate data
# @posts: list of hyper-parameters of the NIG posterior (this corresponds to the output from the function get_posts_my)
# @NNarray: array of nearest neighbours

suppressPackageStartupMessages(library(Matrix))

samp_posts_my_new <- function(posts, NNarray, bayesian = FALSE) {
  # get n, m
  n <- nrow(NNarray)
  m <- ncol(posts[[3]])
  
  # make sure all elements in posts have a place in uhat
  if(ncol(NNarray) < m) {
    stop("The posteriors have more neighbors than the neighbor matrix accounts for!")
  }
  # check if NNarray is an integer matrix; warn if not
  if(! is.integer(NNarray)){
    warning("NNarray should consist of only integers (and/or NAs)!")
  }
  # Chaninging here to output the Xi 's as well
  # make sure posts has 5 elements
  if(length(posts) != 6) {
    stop("There should be 6 elements in the list of posteriors! For simplicity,
         it is recommended to use get_posts to generate these posteriors.")
  }
  # make sure posts elements are of the correct sizes
  if(length(posts[[1]]) != n || length(posts[[2]]) != n || !all(dim(posts[[3]]) == c(n, m)) ||
     !all(dim(posts[[4]]) == c(m, m, n))){
    stop("Please use get_posts to generate the posteriors. Some of 
          the current posteriors have incorrect dimensions.")
  }
  
  if(bayesian){
    n2 <- nrow(NNarray)
    d <- 1/sqrt(rinvgamma(1,posts[[1]][1], posts[[2]][1]))
    
    # New part
    m = min(ncol(posts[[3]]), ncol(NNarray))
    NNarray = NNarray[ , 1:m]
    
    uhat <- Matrix::sparseMatrix(i=1,j=1,x=d, dims=c(n2,n2),triangular=TRUE)
    for(i in 2:n2) {
      gind <- na.omit(NNarray[i,])
      nn <- length(gind)
      d <- rinvgamma(1,posts[[1]][i], posts[[2]][i])
      
      uhat[i,i] <- 1/sqrt(d)
      # posts[[4]][1:nn,1:nn,i] gives U
      # ginv(posts[[4]][1:nn,1:nn,i]) gives U^-1
      # tcrossprod(ginv(posts[[4]][1:nn,1:nn,i])) gives U^-1 (U^-1)^T = (U^T U)^-1 = G_i^-1 in my notation
      uhat[gind,i] <- uhat[i,i]*mvrnorm(1, posts[[3]][i,1:nn], d*tcrossprod(ginv(posts[[4]][1:nn,1:nn,i])))
    }
    return(uhat)
  }else{
    # New part
    m = min(ncol(posts[[3]]), ncol(NNarray))
    NNarray = NNarray[ , 1:m]
    
    # Get posterior mean of d and create the sparse matrix
    d <- (1/sqrt(posts[[2]][1])) * exp(lgamma((2 * posts[[1]][1] + 1)/2) - lgamma(posts[[1]][1]))
    uhat <- sparseMatrix(i = 1, j = 1, x = d, dims = c(n, n), triangular = TRUE)
    for (i in 2:n) {
      # get nearest neighbors and number of them (nn)
      gind <- na.omit(NNarray[i, 1:m])
      nn <- length(gind)
      # Fill in appropriate elements of sparse matrix with posterior means
      uhat[i, i] <- (1/sqrt(posts[[2]][i])) * exp(lgamma((2 * posts[[1]][i] + 1)/2) - lgamma(posts[[1]][i]))
      uhat[gind, i] <- posts[[3]][i, 1:nn] * uhat[i, i]
    }
    return(uhat)
  }
}

# The function to evaluate the log density of matrix-normal distribution MN(0, U, V)
# @X: the value at which the density is to be evaluated (in log-scale)
# @U: Row covariance matrix
# @V: Column covariance matrix

dmatrix_norm <- function(X, U, V){
  n = nrow(X); p = ncol(X)
  
  first = - 0.5 * sum(diag(solve(V, t(X)) %*% solve(U, X)))
  second = -0.5 * n * p * log(2*pi)
  third = - 0.5 * n * determinant(V, logarithm = TRUE)$modulus[1]
  forth = - 0.5 * p * determinant(U, logarithm = TRUE)$modulus[1]
  return(first + second + third + forth)
} 

# The function to suppress messages from inbuilt R Functions (partciularly for adaptMCMC package) 
# @x: the output from an R package
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
