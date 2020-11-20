#' Sieve Estimators
#'
#' For multiple bases - Standard, Bernstein, and Linear Spline
#' Written for data.table arguments.
#' @param Yvar Dependent variable
#' @param Xvar Character vector of control variables
#' @param read_data Of type data.table - source of data
#' @param read_data Of type data.table - points for prediction
#' @param basis Standard for use
#' @param low.bound.Xvar Lower bound of data
#' @param high.bound.Xvar Upper bound of data
#' @param tolerance Tolerance level for matrix inversion
#' @param K degrees of polynomial (or knots for linear spline basis)
#' @keywords sieve polynomial
#' @export
#' @examples
#' my.sieve()

my.sieve <- function(Yvar, Xvar, read_data, predict_data, basis="standard",
                     low.bound.Xvar = 0, high.bound.Xvar = 1,
                     tolerance=1e-16, K=1L) {
  # Select data and add intercept, store obs number for prediction
  .dt <- copy(read_data)[, c(..Yvar, ..Xvar)]
  .dt.pred <- copy(predict_data)
  n_predict <- dim(.dt.pred)[1]

  # Create sieve depending on specified basis
  if (basis=="standard"){
    # Set intercept value and column order, start saving regressor columns
    .dt[,.intercept:=1]
    .dt.pred[,.intercept:=1]
    Xcols <- c(".intercept", Xvar)

    # Create standard basis and save column names
    if (K > 1L) {
      for (d in 2:K) {
        .dt[, paste0(".", Xvar, "_deg_", d) := get(Xvar)^d]
        .dt.pred[, paste0(".", Xvar, "_deg_", d) := get(Xvar)^d]
        Xcols <- c(Xcols, paste0(".", Xvar, "_deg_", d))
      }
    }
  }
  if (basis=="bernstein"){
    # Create Z variable, start saving regressor columns
    .dt[, Z := (get(Xvar)-low.bound.Xvar)/(high.bound.Xvar-low.bound.Xvar)]
    .dt.pred[, Z := (get(Xvar)-low.bound.Xvar)/(high.bound.Xvar-low.bound.Xvar)]
    Xcols <- c()

    # Create Bernstein basis and save column names
    for (k in 0L:K) {
      combination <- factorial(K)/(factorial(K-k)*factorial(k))
      .dt[, paste0("Bz_k_", k) := combination * Z^k * (1-Z)^(K-k)]
      .dt.pred[, paste0("Bz_k_", k) := combination * Z^k * (1-Z)^(K-k)]
      Xcols <- c(Xcols, paste0("Bz_k_", k))
    }
  }
  if (basis=="linear spline"){
    # Set intercept value and column order, start saving regressor columns
    .dt[,.intercept:=1]
    .dt.pred[,.intercept:=1]
    Xcols <- c(".intercept", Xvar)

    # Create truncated power basis and save column names
    for (k in 2L:(K+1)) {
      rkminus1 <- quantile(dt[,..Xvar], probs = (k-1)/(K+1))
      .dt[, paste0("Ind_by_XR_", k, "minus1") :=
           as.numeric(get(Xvar) > rkminus1) * (get(Xvar) - rkminus1)]
      .dt.pred[, paste0("Ind_by_XR_", k, "minus1") :=
           as.numeric(get(Xvar) > rkminus1) * (get(Xvar) - rkminus1)]
      Xcols <- c(Xcols, paste0("Ind_by_XR_", k, "minus1"))
    }
  }

  # Solve OLS and predict conditional mean
  .beta <- (solve(t(.dt[, ..Xcols]) %*%
                    as.matrix(.dt[, ..Xcols]), tol = tolerance)) %*%
      (t(.dt[, ..Xcols]) %*% as.matrix(.dt[, ..Yvar]))
  yhat <- as.matrix(.dt.pred[, ..Xcols]) %*% as.matrix(unlist(.beta))

  return(yhat)
}
