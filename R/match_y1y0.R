#' @title Imputing Missing Potential Outcomes with Matching
#'
#' @description Impute missing potential outcomes for each individual with matching.
#'
#' @details Here are the implementation details for the imputation processes.
#'     Denote \eqn{\hat{Y}^0_i} and \eqn{\hat{Y}^1_i} as the imputed potential
#'     outcomes for individual \eqn{i}. Without loss of generality, if \eqn{A_i = 0}, then
#'     \eqn{\hat{Y}^0_i = Y_i}, and \eqn{\hat{Y}^1_i} is the average of outcomes for the \emph{K} units that are the most
#'     similar to the individual \eqn{i}, i.e.,
#'     \deqn{\hat{Y}_i^0 =  \frac 1 K \sum_{j\in\mathcal{J}_K(i)}Y_j,}
#'     where \eqn{\mathcal{J}_K(i)} represents the set of \eqn{K} matched individuals
#'     with \eqn{A_i = 1}, that are the closest to the individual \eqn{i} in terms of
#'     covariates similarity, and vice versa.
#'
#'
#'
#' @param X A matrix representing covariates, where each row represents the
#'     value of a different covariates for an individual.
#'
#' @param A A vector representing the treatment received by each individual.
#'
#' @param Y A vector representing the observed outcome for each individual.
#'
#' @param K When imputing missing potential outcomes, the average number of
#'     similar individuals are taken based on covariates similarity.
#'
#' @param method The distance measure to be used. It is a argument embed in
#'     \code{\link[stats]{dist}} function.
#'
#' @return Returns a matrix of completed matches, where each row is the imputed \eqn{(Y^1, Y^0)}
#'     for each individual.
#'
#' @export
#'
#' @examples
#' n <- 100
#' p <- 2
#' X <- matrix(rnorm(n*p), ncol = p)
#' A <- sample(c(0,1), n, TRUE)
#' Y <- A * (2*X[,1]) + X[,2]^2 + rnorm(n)
#' match_y1y0(X = X, A = A, Y = Y, K =5)
match_y1y0 <- function(X, A, Y, K = 5, method = 'euclidean'){

  n <- length(A)
  Y1 <- Y0 <- rep(NA, n)
  Y1[A==1] <- Y[A==1]
  Y0[A==0] <- Y[A==0]

  X <- scale(X)
  distance = dist(X, method = method)
  distance = as.matrix(distance)

  sub0 = which(A == 0)
  sub1 = which(A == 1)

  for(j in sub1){
    dd = distance[j,]
    ind1 = sub0[order(dd[sub0])[1:K]]
    Y0[j] <- mean(Y0[ind1])
  }

  for(j in sub0){
    dd = distance[j,]
    ind1 = sub1[order(dd[sub1])[1:K]]
    Y1[j] <- mean(Y1[ind1])
  }

  return(data.frame(Y1 = Y1, Y0 = Y0))
}
