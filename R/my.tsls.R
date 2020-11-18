#' Two-Stage Least Squares
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
#' @keywords tsls linear
#' @export
#' @examples
#' my.tsls()

my.tsls <- function(Yvar, Xvar, Inc, Exc, data, intercept=TRUE,
                   varcov="het", tolerance=1e-16) {
  # Select data and add intercept
  Y <- data[,..Yvar]
  if (Inc != "") {
    X <- data[,c(..Xvar,..Inc)] # Must be char vector
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

  # Coefficients
  P_Z <- as.matrix(Z) %*% solve(t(Z) %*% as.matrix(Z), tol=tolerance) %*% t(Z)
  XtXinv_proj <- solve(t(X) %*% P_Z %*% as.matrix(X), tol=tolerance)
  XtY_proj <- t(X) %*% P_Z %*% as.matrix(Y)
  beta <- XtXinv_proj %*% XtY_proj
  colnames(beta) <- c("coefs")

  # Predicted values, residuals, and squared residuals
  Xhat <- P_Z %*% as.matrix(X)
  Yhat <- as.matrix(X) %*% as.matrix(unlist(beta))
  U <- as.matrix(as.matrix(Y) - Yhat)
  U2 <- c(U) * c(U)

  # Variances (multiple options)
  if (varcov == "het") {
    Sigma <- XtXinv_proj %*% (t(Xhat) %*% diag(U2) %*% as.matrix(Xhat)) %*%
      XtXinv_proj * nrow(X)/(nrow(X)-ncol(X))
  }
  if (varcov == "hom") {
    Sigma <- sum(U2)/(nrow(X)-ncol(X)) * XtXinv_proj
  }

  # Return results
  results <- list(beta = beta, Sigma = Sigma)
  return(results)
}
