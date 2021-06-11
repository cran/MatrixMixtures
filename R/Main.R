#' Fitting for Matrix-Variate Mixture Models
#'
#' Fits, by using expectation-maximization algorithms, mixtures of matrix-variate
#' distributions (normal, t, contaminated normal) to the given data. Can be run
#' in parallel. The Bayesian information criterion (BIC) is used to select the
#' number of groups.
#'
#' @param X A list of dimension \code{N}, where \code{N} is the sample size. Each element of the
#'  list corresponds to an observed p x r matrix.
#' @param G A vector containing the numbers of groups to be tried.
#' @param mod The matrix-variate distribution to be used for the mixture model. Possible
#' values are: \code{"MVN"} for the normal distribution, \code{"MVT"} for the
#' t distribution \code{"MVCN"} for the contaminated normal.
#' @param tol Threshold for Aitken's acceleration procedure. Default value is \code{1.0e-05}.
#' @param maxiter Maximum number of iterations of the algorithms. Default value is \code{10000}.
#' @param ncores A positive integer indicating the number of cores used for running in parallel.
#' Default value is \code{1}.
#' @param verbose Logical indicating whether the running output should be displayed.
#' @return A list with the following elements:
#' \item{flag}{Convergence flag (TRUE - success, FALSE - failure).}
#' \item{pig}{Vector of the estimated mixing proportions (length G).}
#' \item{nu}{Vector of the estimated degree of freedoms (length G). Only for "MVT".}
#' \item{alpha}{Vector of the estimated inliers proportions (length G). Only for "MVCN".}
#' \item{eta}{Vector of the estimated inflation parameters (length G). Only for "MVCN".}
#' \item{M}{Array of the mean matrices (p x r x G).}
#' \item{Sigma}{Array of the estimated row covariance matrices (p x p x G).}
#' \item{Psi}{Array of the estimated column covariance matrices (r x r x G).}
#' \item{class}{Vector of estimated data classification.}
#' \item{z}{Matrix of estimated posterior probabilities (N x G).}
#' \item{v}{Matrix of estimated inlier probabilities (N x G). Only for "MVCN".}
#' \item{lik}{Estimated log-likelihood.}
#' \item{BIC}{Estimated BIC.}
#' @export
#' @importFrom foreach %dopar%
#' @examples
#' data(SimX)
#' res <- MatrixMixt(X = SimX, G = 2, mod = "MVCN")
MatrixMixt <- function(X, G = 1:3, mod, tol = 1.0e-05, maxiter = 10000, ncores = 1, verbose = TRUE) {
  N <- length(X)
  p <- dim(X[[1]])[2]
  n <- dim(X[[1]])[1]
  rd <- N + p * n + 2
  Resall <- list()
  iter <- 1
  l <- NULL
  cluster <- snow::makeCluster(ncores, type = "SOCK")
  doSNOW::registerDoSNOW(cluster)
  for (g in G) {
    if (mod == "MVCN") {
      Res <- foreach::foreach(l = 1:9, .export = c("pig_updateC", "alpha_updateC", "Mg_UpdateC", "Sigma_UpdateC", "Psi_UpdateC", "Eta_UpdateC", "Gauss_Dens_ContC", "Estep_updateC", "EM_ContNorm")) %dopar% {
        tryCatch(EM_ContNorm(X, g, maxiter = maxiter, tol = tol, sm = l * rd, init = "rand"), error = function(e) {
          lt <- list()
          lt[[1]] <- FALSE
          names(lt) <- "flag"
          lt
        })
      }
      Res[[10]] <- tryCatch(EM_ContNorm(X, g, maxiter = maxiter, tol = tol, sm = rd, init = "kmeans"), error = function(e) {
        lt <- list()
        lt[[1]] <- FALSE
        names(lt) <- "flag"
        lt
      })
    } else if (mod == "MVT") {
      Res <- foreach::foreach(l = 1:9, .export = c("pig_updateT", "Mg_UpdateT", "Sigma_UpdateT", "Psi_UpdateT", "nug_updateT", "MVt", "zig_updateT", "Estep2T", "EM_MVT")) %dopar% {
        tryCatch(EM_MVT(X, g, maxiter = maxiter, tol = tol, sm = l * rd, init = "rand"), error = function(e) {
          lt <- list()
          lt[[1]] <- FALSE
          names(lt) <- "flag"
          lt
        })
      }
      Res[[10]] <- tryCatch(EM_MVT(X, g, maxiter = maxiter, tol = tol, sm = rd, init = "kmeans"), error = function(e) {
        lt <- list()
        lt[[1]] <- FALSE
        names(lt) <- "flag"
        lt
      })
    } else if (mod == "MVN") {
      Res <- foreach::foreach(l = 1:9, .export = c("pig_updateN", "Mg_UpdateN", "Sigma_UpdateN", "Psi_UpdateN", "MVN", "zig_updateN", "EM_MVN")) %dopar% {
        tryCatch(EM_MVN(X, g, maxiter = maxiter, tol = tol, sm = l * rd, init = "rand"), error = function(e) {
          lt <- list()
          lt[[1]] <- FALSE
          names(lt) <- "flag"
          lt
        })
      }
      Res[[10]] <- tryCatch(EM_MVN(X, g, maxiter = maxiter, tol = tol, sm = rd, init = "kmeans"), error = function(e) {
        lt <- list()
        lt[[1]] <- FALSE
        names(lt) <- "flag"
        lt
      })
    }
    ConvTest <- lapply(Res, function(x) ConvTestFun(x))
    BICG <- unlist(lapply(ConvTest, function(x) BIC_Calc(x$lik, n, p, N, g, mod)))
    if (max(BICG) == -Inf) {
      Resall[[iter]] <- list(BIC = -Inf)
    } else {
      Restemp <- Res[[which.max(BICG)]]
      Restemp$BIC <- max(BICG)
      Resall[[iter]] <- Restemp
    }
    iter <- iter + 1
  }
  snow::stopCluster(cluster)
  foreach::registerDoSEQ()
  BICall <- unlist(lapply(Resall, function(x) x$BIC))
  if (max(BICall) == -Inf) {
    w.message <- "None of the models could be fitted"
    return(list(flag = F, exit = w.message))
  }
  Gbest <- which.max(BICall)
  if (verbose == TRUE){
  print(paste("The best model chosen by the BIC has ", G[Gbest], " groups", sep = ""))
  }
  if (length(G) > 1) {
    return(Resall[[Gbest]])
  } else {
    return(Resall[[Gbest]])
  }
}
