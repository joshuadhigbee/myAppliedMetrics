#' Local Average Response Function (Abadie, 2003)
#'
#' This function returns point estimates and the variance-covariance matrix.
#' Written for data.table arguments.
#' @param Yvar Dependent variable
#' @param Dvar Treatment variable
#' @param Zvar Instrument
#' @param Xvar1 First-stage running variable (only 1 variable for now)
#' @param Cont1 Controls for first stage (should include FEs or intercept)
#' @param Xvar2 Second stage control variables
#' @param data Of type data.table
#' @param K Degree of polynomial for first stage Xvar1
#' @param tolerance Tolerance level for matrix inversion
#' @keywords ols linear
#' @export
#' @examples
#' my.larf()

my.larf <-
  function(Yvar, Dvar, Zvar, Xvar1, Cont1, Xvar2, data, K=2, tolerance=1e-16) {
  Y <- data[,..Yvar]
  X1 <- data[,c(..Xvar1,..Cont1)] # First stage
  X2 <- data[,..Xvar2] # Second stage
  Z <- data[,..Zvar]
  D <- data[,..Dvar]
  X2 <- X2[,.intercept:=1]
  setcolorder(X2, c(".intercept"))

  # Construct PK (here, written as X1 because it's first stage)
  for (k in 2:K) {
    X1[, paste0("x_", k) := (get(Xvar1))^k]
  }

  # Construct tau_0
  pi_hat <- solve(t(X1) %*% as.matrix(X1), tol=tolerance) %*%
    t(X1) %*% as.matrix(Z)
  tau_0 <- as.matrix(X1) %*% as.matrix(pi_hat)

  # Construct weights
  kappa <- 1 - D*(1-Z)/(1-tau_0) - (1-D)*Z/tau_0
  W <- diag(unlist(kappa))

  # Coefficients
  XtXinv <- solve(t(X2) %*% W %*% as.matrix(X2), tol = tolerance)
  XtY <- t(X2) %*% W %*% as.matrix(Y)
  beta <- XtXinv %*% XtY
  colnames(beta) <- c("coefs")

  # Predicted values, residuals, and squared residuals
  Yhat <- as.matrix(X2) %*% as.matrix(unlist(beta))
  U <- as.matrix(as.matrix(Y) - Yhat)

  # Variances (multiple options)
  N <- nrow(X2)
  dG <- -2 * as.vector(U) * as.matrix(X2)
  M_theta <- matrix(0,nrow=ncol(X2), ncol=ncol(X2))
  for (i in 1:N) {
    k <- unlist(kappa[i])
    M_theta <- M_theta + 1/N * k * 2 * t(X2[i]) %*% as.matrix(X2[i])
  }
  v <- (1-D)*Z/tau_0^2 - D*(1-Z)/(1-tau_0)^2
  delta_X_p1 <- t((as.vector(unlist(v)) * as.matrix(dG))) %*% as.matrix(X1)
  delta_X_p2 <- solve(t(X1) %*% as.matrix(X1), tol=tolerance)
  delta_X <- as.matrix(X1) %*% t(delta_X_p1 %*% delta_X_p2)
  half_meat <- as.vector(unlist(kappa)) * as.matrix(dG) +
    as.matrix(delta_X) * as.vector(unlist(Z - tau_0))
  meat <- t(half_meat) %*% as.matrix(half_meat) / N
  Sigma <- solve(M_theta, tol=tolerance) %*% meat %*%
    solve(M_theta, tol=tolerance) / N

  # Return results
  results <- list(beta = beta, Sigma = Sigma)
  return(results)
}
