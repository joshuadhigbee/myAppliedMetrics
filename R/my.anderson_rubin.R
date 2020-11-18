#' Anderson-Rubin CI for TSLS
#'
#' This function returns point estimates and the variance-covariance matrix.
#' Written for data.table arguments.
#' @param Yvar Dependent variable
#' @param Xvar Endogenous regressor
#' @param Inc Included instruments
#' @param Exc Excluded instruments
#' @param data Of type data.table
#' @param intercept Logical, if including 1 as model intercept
#' @param tolerance Tolerance level for matrix inversion
#' @param beta_low Lower bound for search
#' @param beta_high Upper bound for search
#' @param init_length Size of search grid
#' @param precision Maximum difference between points in and out of CI
#' @keywords anderson-rubin tsls
#' @export
#' @examples
#' my.anderson_rubin()

my.anderson_rubin <-
  function(Yvar, Xvar, Inc, Exc, data, intercept=TRUE,
           tolerance=1e-16, beta_low, beta_high, init_length=10,
           precision) {
  # Select data and add intercept
  Y <- data[,..Yvar]
  if (length(Inc)>0) {
    X <- data[,c(..Xvar,..Inc)] # Must be char vector
    Z <- data[,c(..Inc, ..Exc)] # Must be char vectors
  } else {
    X <- data[,..Xvar] # Must be char vector
    Z <- data[,..Exc] # Must be char vectors
  }
  if (intercept) {
    X <- X[,.intercept:=1]
    Z <- Z[,.intercept:=1]
  }

  # Coefficients
  P_Z <- as.matrix(Z) %*% solve(t(Z) %*% as.matrix(Z), tol=tolerance) %*% t(Z)
  XtXinv_proj <- solve(t(X) %*% P_Z %*% as.matrix(X), tol=tolerance)
  XtY_proj <- t(X) %*% P_Z %*% as.matrix(Y)
  beta_tsls <- XtXinv_proj %*% XtY_proj

  # Useful values to calculate now
  Xhat <- P_Z %*% as.matrix(X)
  ZtZinv <- solve(t(Z) %*% as.matrix(Z), tol=tolerance)

  # Critical value for A-R test, start with null hypothesis for each beta
  crit <- qchisq(0.95, df=ncol(Z))
  beta_grid <-
    data.table(B = seq(from=beta_low, to=beta_high, length.out=init_length),
               reject = 0)

  # Loop through initial beta_grid to find intervals containing lower and
  #   upper bounds for AR CI
  for (b in as.vector(beta_grid[,B])) {
    # Replace endogenous coefficient with hypothesized value b
    beta <- beta_tsls
    beta[1] <- b

    # Predicted values, residuals, and squared residuals
    Yhat <- as.matrix(X) %*% as.matrix(unlist(beta))
    U <- as.matrix(as.matrix(Y) - Yhat)
    U2 <- c(U) * c(U)
    gb <- ZtZinv %*% (t(Z) %*% U)

    # Construct variance, then AR test statistic
    Sigma <- ZtZinv %*% (t(Z) %*% diag(U2) %*% as.matrix(Z)) %*%
      ZtZinv * nrow(Z)/(nrow(Z)-ncol(Z))
    AR <- t(gb) %*% solve(Sigma, tol=tolerance) %*% gb

    # Reject if too large
    if (AR > crit) {
      beta_grid[B==b, reject := 1]
    }
  }

  # Find two grids for searching for lower and upper bounds
    # Lowest and highest accepted beta values
    min_accept <- as.numeric(beta_grid[reject==0, min(B)])
    max_accept <- as.numeric(beta_grid[reject==0, max(B)])
    # Lowest and highest rejected beta values above and below accepted betas
    min_reject <- as.numeric(beta_grid[reject==1 & B > max_accept, min(B)])
    max_reject <- as.numeric(beta_grid[reject==1 & B < min_accept, max(B)])

  # Create search grids for lower and upper bounds
  beta_grid_l <-
    data.table(B = seq(from=max_reject, to=min_accept, length.out=init_length),
               reject = 0)
  beta_grid_u <-
    data.table(B = seq(from=max_accept, to=min_reject, length.out=init_length),
               reject = 0)

  # Loop over each grid to find lower and upper AR CI bounds
  # Lower bound
    diff <- precision*10 # Initialize difference between reject and accept
    while(diff > precision) { # Loop until close enough
    for (b in as.vector(beta_grid_l[,B])) {
      # Replace endogenous coefficient with hypothesized value b
      beta <- beta_tsls
      beta[1] <- b

      # Predicted values, residuals, and squared residuals
      Yhat <- as.matrix(X) %*% as.matrix(unlist(beta))
      U <- as.matrix(as.matrix(Y) - Yhat)
      U2 <- c(U) * c(U)
      gb <- ZtZinv %*% (t(Z) %*% U)

      # Construct variance, then AR test statistic
      Sigma <- ZtZinv %*% (t(Z) %*% diag(U2) %*% as.matrix(Z)) %*%
        ZtZinv * nrow(Z)/(nrow(Z)-ncol(Z))
      AR <- t(gb) %*% solve(Sigma, tol=tolerance) %*% gb

      # Reject if too large
      if (AR > crit) {
        beta_grid_l[B==b, reject := 1]
      }
    }
    # Store "border" values and diff as difference between them
    min_accept_lb <- as.numeric(beta_grid_l[reject==0, min(B)])
    max_reject_lb <- as.numeric(beta_grid_l[reject==1, max(B)])
    diff <- abs(min_accept_lb - max_reject_lb)

    # New grid (if necessary)
    beta_grid_l <- data.table(B = seq(from=max_reject_lb, to=min_accept_lb,
                      length.out=init_length), reject = 0)

    # Print output for monitoring loops
    print(paste("Lower bound iteration.  Difference is", diff))
    }

  # Upper bound
    diff <- precision*10 # Initialize difference between reject and accept
    while(diff > precision) { # Loop until close enough
    for (b in as.vector(beta_grid_u[,B])) {
      # Replace endogenous coefficient with hypothesized value b
      beta <- beta_tsls
      beta[1] <- b

      # Predicted values, residuals, and squared residuals
      Yhat <- as.matrix(X) %*% as.matrix(unlist(beta))
      U <- as.matrix(as.matrix(Y) - Yhat)
      U2 <- c(U) * c(U)
      gb <- ZtZinv %*% (t(Z) %*% U)

      # Construct variance, then AR test statistic
      Sigma <- ZtZinv %*% (t(Z) %*% diag(U2) %*% as.matrix(Z)) %*%
        ZtZinv * nrow(Z)/(nrow(Z)-ncol(Z))
      AR <- t(gb) %*% solve(Sigma, tol=tolerance) %*% gb

      # Reject if too large
      if (AR > crit) {
        beta_grid_u[B==b, reject := 1]
      }
    }
    # Store "border" values and diff as difference between them
    max_accept_ub <- as.numeric(beta_grid_u[reject==0, max(B)])
    min_reject_ub <- as.numeric(beta_grid_u[reject==1, min(B)])
    diff <- abs(max_accept_ub - min_reject_ub)

    # New grid (if necessary)
    beta_grid_u <- data.table(B = seq(from=max_accept_ub, to=min_reject_ub,
                      length.out=init_length), reject = 0)

    # Print output for monitoring loops
    print(paste("Upper bound iteration.  Difference is", diff))
    }

  # Save results
  lower_bound <- min_accept_lb
  upper_bound <- max_accept_ub

  # Return results
  results <- c(lower_bound, upper_bound)
  return(results)
}
