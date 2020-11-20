#' Variance estimator
#'
#' Variance-covariance matrix
#'
#' Written for data.table arguments.
#' @param vars Variables of interest
#' @param data Of type data.table
#' @param tolerance Tolerance level for matrix inversion
#' @keywords variance
#' @export
#' @examples
#' my.variance()

my.variance <- function(vars, data, tolerance=1e-16) {
  X <- data[,..vars] # Xvars must be a character vector with variable names
  n <- dim(X)[1]
  X.demeaned <- X[,lapply(.SD, function(x){x-sum(x)/length(x)})]
  return(t(X.demeaned) %*% as.matrix(X.demeaned)/(n-1))
}
