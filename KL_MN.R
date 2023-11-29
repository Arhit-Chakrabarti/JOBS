KL_MN <- function(data, estimate){
  suppressPackageStartupMessages(library(Matrix))
  # True Sigma, Lambda and number of locations
  Sigma.true = cov2cor(datum$Sigma.true)
  Lambda.true = cov2cor(datum$Lambda.true)
  n_locs = nrow(datum$locs)
  # Posterior estimates of Lambda and Sigma
  Lambda_post = cov2cor(estimate$Lambda_post)
  Sigma_post = cov2cor(solve(t(estimate$U_post),  solve(estimate$U_post)))
  num_reps = nrow(Lambda.true)
  # Calculate the KL divergence for Matrix Normal 
  V1 = Sigma.true; V2 = Sigma_post; U1 = Lambda.true; U2 = Lambda_post  # Look at notes for definitions of U1, U2, V1, V2
  # KL divergence calculation. All log(det) are calculated using determinant( , log = TRUE) and appropriate ratios are replaced by differences. Formula is complicated. Be sure to double check
  KL_Matrix_Normal = log((sum(diag(solve(V2, V1))) * sum(diag(solve(U2, U1))) - (num_reps*(determinant(V1, log = TRUE)$modulus - determinant(V2, log = TRUE)$modulus))[1] - (n_locs * (determinant(U1, log = TRUE)$modulus - determinant(U2, log = TRUE)$modulus))[1] - (n_locs * num_reps))/2)
  return(KL_Matrix_Normal)
}