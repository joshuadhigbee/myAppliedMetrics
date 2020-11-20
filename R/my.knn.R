#' K-nearest neighbors smoother
#'
#' For multiple metrics - Euclidean and Mahalanobis
#' Written for data.table arguments.
#' @param Yvar Dependent variable
#' @param Xvars Character vector of control variables
#' @param read_data Of type data.table - source of data
#' @param predict_data Of type data.table - points for prediction
#' @param metric Metric for comparison
#' @param tolerance Tolerance level for matrix inversion
#' @param K number of neighbors
#' @keywords nearest neighbors
#' @export
#' @examples
#' my.knn()

my.knn <- function(Yvar, Xvars, read_data, predict_data, metric="mahalanobis",
                     tolerance=1e-16, K=1L) {
  # Select data and add intercept, store obs number
  .dt <- copy(read_data)[, c(..Yvar, ..Xvars)]
  .dt.pred <- copy(predict_data)
  n_predict <- dim(.dt.pred)[1]

  # Create var-cov matrix for distance using specified metrics
  if (metric=="mahalanobis"){
    Sigma_inv <- solve(.dt[, my.variance(..Xvars, .SD, tolerance = tolerance)],
                       tolerance = tolerance)
  }
  if (metric=="euclidean"){
    Sigma_inv <- diag(length(Xvars))
  }

  # Create empty vector for storing conditional mean estimates
  yhat <- rep(1, n_predict)

  # Match each observation and take mean among K nearest neighbors
  for (i_row in 1:n_predict) {
    # Create dataset of other observations and distance from observation
    r <- unlist(.dt.pred[i_row, ..Xvars])
    mat <- as.matrix(sweep(.dt[, ..Xvars], 2, r))
    .dt[, .distance := sqrt(rowSums( (mat %*% Sigma_inv) * mat))]

    # Find K smallest distances after ordering by distance
    setorder(.dt, .distance)
    yhat[i_row] <- unlist(.dt[1L:K, my.mean(..Yvar, .SD)])
  }

  return(yhat)
}
