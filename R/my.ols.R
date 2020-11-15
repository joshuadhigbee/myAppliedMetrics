#' Ordinary Least Squares
#'
#' This function returns point estimates and the variance-covariance matrix.
#' Written for data.table arguments.
#' @param Yvar Dependent variable
#' @param Xvars Character vector of control variables
#' @param data Of type data.table
#' @param intercept Logical, if including 1 as model intercept
#' @param varcov Method for variance-covariance estimation
#' @param tolerance Tolerance level for matrix inversion
#' @keywords ols linear
#' @export
#' @examples
#' my.ols()

my.ols <- function(Yvar, Xvars, data, intercept=TRUE,
                   varcov="het", tolerance=1e-16) {
  # Select data and add intercept
  Y <- data[,..Yvar]
  X <- data[,..Xvars] # Must be char vector
  if (intercept) {
    X <- X[,.intercept:=1]
    setcolorder(X, c(".intercept"))
  }

  # Coefficients
  XtXinv <- solve(t(X) %*% as.matrix(X), tol = tolerance)
  XtY <- t(X) %*% as.matrix(Y)
  beta <- XtXinv %*% XtY
  colnames(beta) <- c("coefs")

  # Predicted values, residuals, and squared residuals
  Yhat <- as.matrix(X) %*% as.matrix(unlist(beta))
  U <- as.matrix(as.matrix(Y) - Yhat)
  U2 <- c(U) * c(U)

  # Variances (multiple options)
  if (varcov == "het") {
    Sigma <- XtXinv %*% (t(X) %*% diag(U2) %*% as.matrix(X)) %*% XtXinv *
      nrow(X)/(nrow(X)-ncol(X))
  }
  if (varcov == "hom") {
    Sigma <- sum(U2)/(nrow(X)-ncol(X)) * XtXinv
  }

  # Return results
  results <- list(beta = beta, Sigma = Sigma)
  return(results)
}
