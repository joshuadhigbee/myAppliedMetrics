#' Jackknife IV Estimator
#'
#' This function returns point estimates and the variance-covariance matrix.
#' Written for data.table arguments.
#' @param Yvar Dependent variable
#' @param Xvar Endogenous variable
#' @param Inc Included instruments (character vector of controls)
#' @param Exc Excluded instruments (character vector)
#' @param data Of type data.table
#' @param intercept Logical, if including 1 as model intercept
#' @param varcov Method for variance-covariance estimation
#' @param tolerance Tolerance level for matrix inversion
#' @keywords ols linear
#' @export
#' @examples
#' my.jackknife_iv()

my.jackknife_iv <- function(Yvar, Xvar, Inc, Exc, data, intercept=TRUE,
                   varcov="het", tolerance=1e-16) {
  # Select data and add intercept
  Y <- data[,..Yvar]
  if (Inc!="") { # If there are included instruments besides intercept
    X <- data[,c(..Inc, ..Xvar)] # Must be char vector
    Z <- data[,c(..Inc, ..Exc)] # Must be char vectors
  } else {
    X <- data[,..Xvar] # Must be char vector
    Z <- data[,..Exc] # Must be char vectors
  }
  if (intercept) {
    X <- X[,.intercept:=1]
    setcolorder(X, c(".intercept"))
    Z <- Z[,.intercept:=1]
    setcolorder(Z, c(".intercept"))
  }

  # Construct jackknife instrument
  ZtZinv <- solve(t(Z) %*% as.matrix(Z), tol=tolerance)
  h_i <- diag(as.matrix(Z) %*% ZtZinv %*% t(Z))
  Pi_hat <- t(X) %*% as.matrix(Z) %*% ZtZinv
  X_jack <- (1-h_i)^(-1)*(as.matrix(Z) %*% t(Pi_hat)  - h_i * X)

  # Coefficients
  XtXinv_jack <- solve(t(X_jack) %*% as.matrix(X), tol=tolerance)
  XtY_jack <- t(X_jack) %*% as.matrix(Y)
  beta <- XtXinv_jack %*% XtY_jack
  colnames(beta) <- c("coefs")

  # Predicted values, residuals, and squared residuals
  Yhat <- as.matrix(X) %*% as.matrix(unlist(beta))
  U <- as.matrix(as.matrix(Y) - Yhat)
  U2 <- c(U) * c(U)

  # Variances (multiple options)
  if (varcov == "het") {
    Sigma <- XtXinv_jack %*% (t(X_jack) %*% diag(U2) %*% as.matrix(X_jack)) %*%
      t(XtXinv_jack) * nrow(X)/(nrow(X)-ncol(X))
  }
  if (varcov == "hom") {
    Sigma <- sum(U2)/(nrow(X)-ncol(X)) * XtXinv_jack %*%
      (t(X_jack) %*% as.matrix(X_jack)) %*% t(XtXinv_jack)
  }

  # Return results
  results <- list(beta = beta, Sigma = Sigma)
  return(results)
}
