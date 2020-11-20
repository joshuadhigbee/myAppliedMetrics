#' K-nearest neighbors matching - binary treatment
#'
#' For multiple metrics - Euclidean and Mahalanobis
#' Written for data.table arguments.
#' @param Yvar Dependent variable
#' @param Xvars Character vector of control variables
#' @param binary_treat_var Binary treatment indicator
#' @param read_data Of type data.table - source of data
#' @param predict_data Of type data.table - points for prediction
#' @param metric Metric for comparison
#' @param tolerance Tolerance level for matrix inversion
#' @param K number of neighbors
#' @keywords nearest neighbors
#' @export
#' @examples
#' my.knn.binary()

my.knn.binary <- function(Yvar, Xvars, read_data,
          metric="mahalanobis", binary_treat_var, tolerance=1e-16, K=1L) {
  # Select data and add intercept, store obs number
  .dt <- copy(read_data)[, c(..Yvar, ..Xvars, ..binary_treat_var)]
  n_predict <- dim(.dt)[1]

  # Create var-cov matrix for distance using specified metrics
  if (metric=="mahalanobis"){
    Sigma_inv <- solve(.dt[get(binary_treat_var)==0,
                           my.variance(..Xvars, .SD, tolerance = tolerance)],
                       tolerance = tolerance)
  }
  if (metric=="euclidean"){
    Sigma_inv <- diag(length(Xvars))
  }

  # Create vector for storing conditional mean estimates
  y1_hat <- rep(1, n_predict)
  y0_hat <- rep(1, n_predict)

  # Create dataset for number of times matched and index for rows in data
  Km <- data.table(times_matched = rep(0, n_predict),
                   index = 1:n_predict)
  .dt[, index := 1:n_predict]

  for (i_row in 1:n_predict) {
    # Check treatment status of each row, look at opposite-value treatment obs
    r_treat_status <- unlist(.dt[i_row, ..binary_treat_var])
    .dt.temp <- copy(.dt)[get(binary_treat_var)==1-r_treat_status]

    # Create dataset of other observations and distance from observation
    r <- unlist(.dt[i_row, ..Xvars])
    mat <- as.matrix(sweep(.dt.temp[, ..Xvars], 2, r))
    .dt.temp[, .distance := sqrt(rowSums( (mat %*% Sigma_inv) * mat))]

    # Find K smallest distances after ordering by distance (keep index)
    setorder(.dt.temp, .distance)
    yhat_others <- unlist(.dt.temp[1L:K, my.mean(..Yvar, .SD)])
    matched_obs_index <- unlist(.dt.temp[1L:K, index])

    # Predict by treatment status
    if (r_treat_status == 0) {
      y0_hat[i_row] <- .dt[i_row, ..Yvar]
      y1_hat[i_row] <- yhat_others
    }
    if (r_treat_status == 1) {
      y0_hat[i_row] <- yhat_others
      y1_hat[i_row] <- .dt[i_row, ..Yvar]
    }

    # Store matched obs
    for (ind in matched_obs_index) {
      Km[index==ind, times_matched := times_matched + 1]
    }
  }

  # Estimate parameters of interest
  .dt[, y1_hat_est := unlist(y1_hat)]
  .dt[, y0_hat_est := unlist(y0_hat)]
  ATE <- .dt[, ATE_i := y1_hat_est - y0_hat_est][, my.mean(c("ATE_i"), .SD)]
  ATT <- .dt[, ATT_i := y1_hat_est - y0_hat_est
             ][get(binary_treat_var)==1, my.mean(c("ATT_i"), .SD)]

  # Empty sigma(X) vector
  sigmaX <- rep(0, n_predict)

  # Variance estimation
  for (i_row in 1:n_predict) {
    # Check treatment status of each row, look at opposite-value treatment obs
    r_treat_status <- unlist(.dt[i_row, ..binary_treat_var])
    .dt.temp <- copy(.dt)[get(binary_treat_var)==r_treat_status]

    # Create dataset of other observations and distance from observation
    r <- unlist(.dt[i_row, ..Xvars])
    mat <- as.matrix(sweep(.dt.temp[, ..Xvars], 2, r))
    .dt.temp[, .distance := sqrt(rowSums( (mat %*% Sigma_inv) * mat))]

    # Find K smallest distances after ordering by distance
    setorder(.dt.temp, .distance)
    yhat_others <- unlist(.dt.temp[1L:K, my.mean(..Yvar, .SD)])

    # Estimate Ybar
    Ybar <- 1/(K+1)*(unlist(.dt[i_row, ..Yvar]) + K*yhat_others)

    # Estimate sigma(X_i)
    sum <- (unlist(.dt[i_row, ..Yvar])-Ybar)^2
    for (m in 1:K) {
      sum <- sum + (unlist(.dt.temp[1L:K, ..Yvar][m])-Ybar)^2
    }

    sigmaX[i_row] <- sum/K
  }

  # Return variance using formulas in Abadie et al.
  Km[, sigmaX := unlist(sigmaX)]
  Km[, W := unlist(.dt[, ..binary_treat_var])]
  Km[, K_mi := times_matched/(K)]
  Km[, K_mi2 := times_matched/(K)^2]
  Km[, sum_term_ATE := ((1+K_mi)^2)*sigmaX] # Replace times_matched with K_mi
  Km[, sum_term_ATT := ((W - (1-W)*K_mi)^2)*sigmaX]
  Var_ATE <- (1/.dt[,.N]^2)*sum(Km[,sum_term_ATE])
  Var_ATT <- (1/.dt[get(binary_treat_var)==1,.N]^2)*sum(Km[,sum_term_ATT])

  # Return results
  results <- list(ATE = ATE, ATT = ATT, Var_ATE = Var_ATE, Var_ATT = Var_ATT,
                  N_obs = n_predict)
  return(results)
}
