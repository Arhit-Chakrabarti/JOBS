suppressPackageStartupMessages(source("gen_data.R"))
# Function to run our MCMC algorithm
# @data: a list consisting of the spatial data, spatial locations, true log-likelihood etc. The data is from the generate_data function
# @m: the number of neighbours
# @n_iterations: number of MCMC iterations.
# @burn: number of samples to be considered as burn-in
# @thin: thinning factor of the MCMC samples
posterior_est <- function(data, m = 10, n_iterations = 1000, burn = 500, thin = 1){
  if(burn >= n_iterations){
    stop("Number of iterations is less than number of burn-in")
  }
  dataa <- data$data
  locs <- data$locs
  
  # Parameters of the data
  num_reps <- nrow(dataa)
  n_locs <- ncol(dataa)
  
  ##### PRE-PROCESSING THE DATA
  # Compute maximin ordering
  order <- order_maximin_dist(fields::rdist(locs))
  
  # Reorder data and location by maximin ordering
  dataa <- dataa[, order]
  locs <- locs[order, ]
  # Find the Euclidean neighbors
  nearest_neighbors <- find_nn_dist(fields::rdist(locs), n_locs)
  
  N = nrow(dataa) # Define the number of replicates

  
  nu = ncol(t(dataa)) # Degrees of freedom for the IW prior
  Psi = diag(1, num_reps) # Set scaling parameter for IW prior
  Lambda = as.matrix(forceSymmetric(riwish(v = num_reps, S = Psi))) # Starting value for the IW sampling (taken to be equal to true value)
  
  # Define lists to store the samples
  U_samples <- list()
  d_samples <- list()
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
    thetas_new <- quiet(adaptMCMC::MCMC(minus_loglikeli_my_new, datum = dataa, 
                                        NNarray = nearest_neighbors,
                                        Lambda = Lambda,
                                        m = m,
                                        init = init_theta, negativ = FALSE,
                                        scale = scale_mat, adapt = TRUE, 
                                        acc.rate = 0.234, n = nruns, 
                                        showProgressBar = FALSE)$samples[nruns, ])
    
    theta_samples[i, ] = thetas_new
    # Convert the hyper-parameters thetas to priors which is a function of Lambda from previous iteration
    priors = thetas_to_priors_my(thetas_new, n = nrow(nearest_neighbors), m = m)
    # get initial posterior sample
    post = get_posts_my(datum = dataa, priors = priors, Lambda = Lambda, NNarray = nearest_neighbors)
    
    # Sample posterior sparse matrix U as a function of Lambda from previous iteration
    U_samples[[i]] = samp_posts_my_new(post, NNarray = nearest_neighbors, bayesian = TRUE)
    # Store the diagonal elements of U
    d_samples[[i]] = diag(U_samples[[i]])
    
    SS = matrix(0, nrow = N, ncol = N) # Define the sample sum of square matrix
    # ALTERNATE
    for(ind in 1:n_locs){
      if(ind == 1){
        mean.normal = rep(0, num_reps)
        d = 1/((d_samples[[i]][ind])^2)
        SS = SS + tcrossprod(dataa[, ind])/d
      }else{
        X = as.matrix(post[[5]][[ind]])
        d = 1/((d_samples[[i]][ind])^2)
        u =  sqrt(d) * U_samples[[i]][as.numeric(post[[6]][[ind]]), ind]  
        mean.normal = as.numeric(X %*% u)
        SS = SS + (tcrossprod(dataa[, ind] - mean.normal))/d
      }
    }
    # Draw sample from IW i.e. from full conditional [Lambda|-]
    Lambda = riwish(v = nu + n_locs, S = as.matrix(forceSymmetric(Psi + SS)))
    # Store Lambda as sample
    Lambda_samples[[i]] = as.matrix(forceSymmetric(Lambda))
    # current theta becomes old theta for next itertaion
    # thetas_old = thetas_new
    log_like <- 0
    for(ind in 1:n_locs){
      if(ind == 1){
        mean.normal = rep(0, num_reps)
        d = 1/((d_samples[[i]][ind])^2)
      }else{
        X = as.matrix(post[[5]][[ind]])
        d = 1/((d_samples[[i]][ind])^2)
        u =  sqrt(d) * U_samples[[i]][as.numeric(post[[6]][[ind]]), ind]  
        mean.normal = as.numeric(X %*% u)
        
      }
      
      log_like[ind] = dmvnorm(x = dataa[, ind], mean = mean.normal, sigma = d * Lambda_samples[[i]], log = TRUE)
    }
    log_likelihood[[i]] = sum(log_like)
  }# End of Gibbs Sampling
  end_time <- Sys.time()
  time.taken <- end_time - start_time
  
  U_post = Matrix(0, nrow = n_locs, ncol = n_locs, sparse = TRUE)
  Lambda_post = matrix(0, nrow = num_reps, ncol = num_reps)
  theta_post = 0

  samples <- seq(from = (burn + 1), to = n_iterations, by = thin)
  
  for(i in samples){
    U_post = U_post + U_samples[[i]]  # Sum of all matrices from Posterior of U
    Lambda_post = Lambda_post + Lambda_samples[[i]] # Sum of all matrices from Posterior of Lambda
  }
  
  U_post = U_post/length(samples) # Posterior Mean as an estimate of U
  Lambda_post = Lambda_post/length(samples) # Posterior Mean as an estimate of Lambda
  theta_post = apply(theta_samples, MARGIN = 2, FUN = mean)
  
  log_like_new <- list()
  
  for(i in (burn + 1):n_iterations){
    log_like_new[[i]] <- dmatrix_norm(X = dataa, U = cov2cor(Lambda_samples[[i]]), V = cov2cor(solve(t(U_samples[[i]]),  solve(U_samples[[i]]))))
  }
  
  ll <- data.frame(Iterations = 1:length(samples), LL = unlist(log_likelihood[samples]))
  suppressPackageStartupMessages(library(tidyverse))
  LL.plot <- ll %>% ggplot(aes(x = Iterations, y = LL)) + geom_line() +
    labs(title = "Log-likelihood", x = "Iterations post burn-in", y = "") + 
    theme(axis.text=element_text(size=12),                                                                          axis.title=element_text(size=12), 
          plot.title = element_text(size=12)) + geom_abline(intercept = data$ll.true, slope = 0, col = "red")

  
  suppressPackageStartupMessages(library(forecast))
  ACF.plot <-  ggAcf(x = unlist(log_likelihood[samples]), lag.max = 40) + ggtitle("ACF of log-likelihood") + labs(y = "") + ylim(c(-0.35, 0.65)) + theme(axis.text=element_text(size=12),                                                                           axis.title=element_text(size=12))
                                                                                                                                  
 return(list(U_post = U_post,
             Lambda_post = Lambda_post,
             LL.plot = LL.plot,
             ACF.plot = ACF.plot,
             log_like_new = unlist(log_like_new[samples])))
}





