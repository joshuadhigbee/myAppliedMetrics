#' Local polynomial estimator
#'
#' This function returns point estimates and the variance-covariance matrix.
#' Written for data.table arguments.
#' @param Yvar Dependent variable
#' @param Xvar Running variable
#' @param read_data Of type data.table - data for reading
#' @param predict_data Of type data.table - Xvar points for predicting values
#' @param degree Degree of polynomial
#' @param bandwidth Bandwidth of data for Xvar
#' @param tolerance Tolerance level for matrix inversion
#' @keywords local linear
#' @export
#' @examples
#' my.local_poly()

my.local_poly <- function(Yvar, Xvar, read_data, predict_data,
                          tolerance=1e-16, degree=1L, bandwidth=1) {
  # Select data and add intercept, store obs number for prediction
  .dt <- copy(read_data)[,.intercept:=1]
  .dt.pred <- copy(predict_data)[,.intercept:=1]
  n_predict <- dim(.dt.pred)[1]

  # Create terms for polynomial degree N > 1 (default = 1)
  Xcols <- c(".intercept") # Keep track of X variable names
  if (degree>=1L) { Xcols <- c(Xcols, Xvar)}
  if (degree > 1L) {
    for (d in 2L:degree) {
      .dt[, paste0(Xvar, "_", d) := get(Xvar)^d]
      .dt.pred[, paste0(Xvar, "_", d) := get(Xvar)^d]
      Xcols <- c(Xcols, paste0(Xvar, "_", d))
    }
  }

  # Create vector for storing conditional mean estimates
  yhat <- rep(1, n_predict)

  # Loop through observations to create mean at each prediction point
  for (i_row in 1:n_predict) {
    # Prediction point
    x_point <- unlist(.dt.pred[i_row, ..Xvar])

    # Subset data for local estimation
    .dt.temp <- .dt[get(Xvar) %between%
                    c(x_point - bandwidth, x_point + bandwidth)]
    .Y <- .dt.temp[,..Yvar]
    .X <- .dt.temp[,..Xcols]

    # Solve OLS and predict conditional mean
    .beta <- (solve(t(.X) %*% as.matrix(.X), tol = tolerance)) %*%
      (t(.X) %*% as.matrix(.Y))
    yhat[i_row] <- as.matrix(.dt.pred[,..Xcols][i_row]) %*%
      as.matrix(unlist(.beta))
  }

  return(yhat)
}
