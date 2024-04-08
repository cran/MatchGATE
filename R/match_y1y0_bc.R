#' @title Imputing Missing Potential Outcomes with Bias-Corrected Matching
#'
#' @description Impute missing potential outcomes for each individual with bias-corrected matching.
#'
#' @details Here are the implementation details for the imputation processes.
#'     Denote \eqn{\hat{Y}^0_i} and \eqn{\hat{Y}^1_i} as the imputed potential
#'     outcomes for individual \eqn{i}. For example, if \eqn{A_i = 0}, then \eqn{\hat{Y}^0_i =  Y^0_i}.
#'     However, for obtaining \eqn{\hat{Y}^1_i}, we require to introduce an outcome
#'     regression function \eqn{\mu_1(X)} for \eqn{Y^1}. Let \eqn{\hat{\mu}_1(X)} be the fitted value of
#'     \eqn{\mu_1(X)}, then \eqn{\hat{Y}^1_i} is defined as follows,
#'     \deqn{\hat{Y}_i^1 =  \frac 1 K \sum_{j\in\mathcal{J}_K(i)}\{Y_j+
#'     \hat{\mu}_1(X_i)-\hat{\mu}_1(X_j)\},}
#'     where \eqn{\mathcal{J}_K(i)} represents the set of \eqn{K} matched individuals
#'     with \eqn{A_i = 1}, that are the closest to the individual \eqn{i} in terms of
#'     covariates similarity, and vice versa.
#'
#' @param X A matrix representing covariates, where each row represents the
#'     value of a different covariates for an individual.
#'
#' @param A A vector representing the treatment received by each individual.
#'
#' @param Y A vector representing the observed outcome for each individual.
#'
#'
#' @param miu1.hat The estimated outcome regression function for \eqn{Y^1}.
#'
#' @param miu0.hat The estimated outcome regression function for \eqn{Y^0}.
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
#' n = 100
#' X1 <- runif(n, -0.5,0.5)
#' X2 <- sample(c(0,1,2), n, TRUE)
#' X = cbind(X1, X2)
#' A = sample(c(0,1), n, TRUE)
#' Y = A * (2*X1) + X1 + X2^2 + rnorm(n)
#' miu1_hat <- cbind(1,X) %*% as.matrix(lm(Y ~ X, subset = A==1)$coef)
#' miu0_hat <- cbind(1,X) %*% as.matrix(lm(Y ~ X, subset = A==0)$coef)
#' match_y1y0_bc(X = X, A = A, Y = Y, miu1.hat = miu1_hat,
#'               miu0.hat = miu0_hat, K = 5)
#'
#'
match_y1y0_bc <- function(X, A, Y, miu1.hat, miu0.hat,
                          K = 5, method = 'euclidean'){

  n <- length(A)
  Y1 <- Y0 <- rep(NA, n)
  Y1[A==1] <- Y[A==1]
  Y0[A==0] <- Y[A==0]

  X <- scale(X)
  distance = dist(X, method = 'euclidean')
  distance = as.matrix(distance)

  sub0 = which(A == 0)
  sub1 = which(A == 1)

  for(j in sub1){
    dd = distance[j,]
    ind1 = sub0[order(dd[sub0])[1:K]]
    Y0[j] <- miu0.hat[j] + mean(Y0[ind1] - miu0.hat[ind1])
  }

  for(j in sub0){
    dd = distance[j,]
    ind1 = sub1[order(dd[sub1])[1:K]]
    Y1[j] <- miu1.hat[j] + mean(Y1[ind1] - miu1.hat[ind1])
  }

  return(data.frame(Y1 = Y1, Y0 = Y0))
}
