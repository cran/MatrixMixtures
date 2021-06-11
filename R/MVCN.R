pig_updateC <- function(z, N) {
  pig <- colSums(z) / N
  return(pig)
}

alpha_updateC <- function(z, v, N) {
  G <- ncol(z)
  alpha <- NULL
  for (g in 1:G) {
    alpha[g] <- max(0.51, sum(z[, g] * v[, g]) / sum(z[, g]))
  }
  return(alpha)
}

Mg_UpdateC <- function(X, z, v, eta) {
  Mg <- list()
  G <- ncol(z)
  N <- length(X)
  for (g in 1:G) {
    stemp <- (z * v)[, g] + (z * (1 - v))[, g] / eta[g]
    numer <- Reduce("+", lapply(1:N, function(i) stemp[i] * X[[i]]))
    Mg[[g]] <- numer / sum(stemp)
  }
  return(Mg)
}

Sigma_UpdateC <- function(X, z, v, eta, M, Psiin) {
  Sigmain <- list()
  p <- dim(Psiin[[1]])[1]
  detSig <- NULL
  G <- ncol(z)
  N <- nrow(z)
  for (g in 1:G) {
    stemp <- (z * v)[, g] + (z * (1 - v))[, g] / eta[g]
    numer <- Reduce("+", lapply(1:N, function(i) stemp[i] * (X[[i]] - M[[g]]) %*% Psiin[[g]] %*% t(X[[i]] - M[[g]])))
    temp <- numer / (p * colSums(z)[g])
    Sigmain[[g]] <- solve(temp)
    detSig[g] <- det(Sigmain[[g]])
  }
  return(list(Sigmain, detSig, prob = 0))
}

Psi_UpdateC <- function(X, z, v, eta, M, Sigmain) {
  Psiin <- list()
  n <- dim(Sigmain[[1]])[1]
  detPsi <- NULL
  G <- ncol(z)
  N <- nrow(z)
  for (g in 1:G) {
    stemp <- (z * v)[, g] + (z * (1 - v))[, g] / eta[g]
    numer <- Reduce("+", lapply(1:N, function(i) stemp[i] * t(X[[i]] - M[[g]]) %*% Sigmain[[g]] %*% (X[[i]] - M[[g]])))
    temp <- numer / (n * colSums(z)[g])
    Psiin[[g]] <- solve(temp)
    detPsi[g] <- det(Psiin[[g]])
  }
  return(list(Psiin, detPsi, prob = 0))
}

Eta_UpdateC <- function(X, M, Sigmain, Psiin, z, v) {
  n <- dim(Sigmain[[1]])[1]
  p <- dim(Psiin[[1]])[1]
  G <- length(M)
  N <- length(X)
  etanew <- NULL
  zv <- z * (1 - v)
  for (g in 1:G) {
    denom <- n * p * sum(zv[, g])
    etatemp <- Reduce("+", lapply(1:N, function(i) zv[i, g] * sum(diag(Sigmain[[g]] %*% (X[[i]] - M[[g]]) %*% Psiin[[g]] %*% t(X[[i]] - M[[g]]))))) / denom
    etanew[g] <- max(c(1.001, etatemp))
  }
  return(etanew)
}

Gauss_Dens_ContC <- function(X, M, Sigmain, Psiin, detSig, detPsi, eta, alpha) {
  n <- dim(Sigmain)[1]
  p <- dim(Psiin)[1]
  logdensNorm <- -(n * p / 2) * log(2 * pi) + (p / 2) * log(detSig) + (n / 2) * log(detPsi) - 0.5 * sum(diag(Sigmain %*% (X - M) %*% Psiin %*% t(X - M)))
  logdensCont <- -(n * p / 2) * log(2 * pi) + (p / 2) * log(detSig) + (n / 2) * log(detPsi) - (1 / (2 * eta)) * sum(diag(Sigmain %*% (X - M) %*% Psiin %*% t(X - M))) - (n * p / 2) * log(eta)
  densNorm <- alpha * exp(logdensNorm)
  densCont <- (1 - alpha) * exp(logdensCont)
  Density <- densNorm + densCont
  return(list(Density, densNorm, densCont))
}

Estep_updateC <- function(X, pig, alpha, eta, M, Sigmain, Psiin, detSig, detPsi) {
  G <- length(M)
  N <- length(X)
  zmat <- matrix(0, nrow = N, ncol = G)
  vmat <- matrix(0, nrow = N, ncol = G)
  for (g in 1:G) {
    temp <- lapply(X, function(x) Gauss_Dens_ContC(x, M[[g]], Sigmain[[g]], Psiin[[g]], detSig[g], detPsi[g], eta[g], alpha[g]))
    zmat[, g] <- unlist(lapply(temp, function(x) x[[1]] * pig[g]))
    vmat[, g] <- unlist(lapply(temp, function(x) x[[2]] / x[[1]]))
  }
  return(list(zmat, vmat))
}

EM_ContNorm <- function(X, G, maxiter = 1000, tol = 1.0e-03, sm = 1, init = "rand") {
  n <- dim(X[[1]])[1]
  p <- dim(X[[1]])[2]
  N <- length(X)
  M <- list()
  Sigmain <- list()
  Psiin <- list()
  detSig <- NULL
  detPsi <- NULL
  if (init == "kmeans") {
    vecX <- matrix(Reduce("rbind", lapply(X, function(x) as.vector(x))), N, n * p)
    zinit <- matrix(0, N, G)
    withr::with_seed(sm, kmeansX <- stats::kmeans(vecX, centers = G, nstart = 25)$cluster)
    for (i in 1:N) {
      zinit[i, kmeansX[i]] <- 1
    }
  } else if (init == "rand") {
    withr::with_seed(sm, zinit <- matrix(stats::runif(N * G), nrow = N))
    zinit <- zinit / rowSums(zinit)
  }
  z <- zinit
  v <- matrix(0.99, nrow = N, ncol = G)
  for (g in 1:G) {
    M[[g]] <- Reduce("+", lapply(1:N, function(i) z[i, g] * X[[i]])) / sum(z[, g])
    Sigma <- Reduce("+", lapply(1:N, function(i) z[i, g] * (X[[i]] - M[[g]]) %*% t(X[[i]] - M[[g]]))) / (p * sum(z[, g]))
    Sigmain[[g]] <- solve(Sigma)
    detSig[g] <- det(Sigmain[[g]])
    Psi <- Reduce("+", lapply(1:N, function(i) z[i, g] * t(X[[i]] - M[[g]]) %*% (X[[i]] - M[[g]]))) / (p * sum(z[, g]))
    Psiin[[g]] <- solve(Psi)
    detPsi[g] <- det(Psiin[[g]])
    eta <- rep(2, G)
  }
  conv <- 0
  prob <- 0
  iter <- 1
  lik <- NULL
  while (conv == 0) {
    pig <- pig_updateC(z, N)
    alpha <- alpha_updateC(z, v, N)
    M <- Mg_UpdateC(X, z, v, eta)
    eta <- Eta_UpdateC(X, M, Sigmain, Psiin, z, v)
    Sigmatest <- Sigma_UpdateC(X, z, v, eta, M, Psiin)
    if (Sigmatest$prob == 1) {
      return(list(flag = F))
    }
    Sigmain <- Sigmatest[[1]]
    detSig <- Sigmatest[[2]]
    Psitest <- Psi_UpdateC(X, z, v, eta, M, Sigmain)
    if (Psitest$prob == 1) {
      return(list(flag = F))
    }
    Psiin <- Psitest[[1]]
    detPsi <- Psitest[[2]]
    Estep <- Estep_updateC(X, pig, alpha, eta, M, Sigmain, Psiin, detSig, detPsi)
    zmat <- Estep[[1]]
    z <- zmat / rowSums(zmat)
    v <- Estep[[2]]
    if (any(is.nan(v))) {
      return(list(flag = F))
    }
    lik[iter] <- sum(log(rowSums(zmat)))
    if (is.na(lik[iter])) {
      w.message <- "NA in likelihood!"
      return(list(flag = F, exit = w.message))
    }
    if (is.infinite(lik[iter])) {
      w.message <- "Infinite likelihood!"
      return(list(flag = F, exit = w.message))
    }
    if (iter > 1) {
      if ((lik[iter] - lik[iter - 1]) < 0) {
        w.message <- paste("Decreasing Likelihood! G=", G, " after ", iter, " iterations", sep = "")
        return(list(flag = F, exit = w.message))
      }
    }
    if (iter > 5) {
      if ((lik[iter - 1] - lik[iter - 2]) == 0) {
        conv <- 1
      } else {
        ak <- (lik[iter] - lik[iter - 1]) / (lik[iter - 1] - lik[iter - 2])
        linf <- lik[iter - 1] + (lik[iter] - lik[iter - 1]) / (1 - ak)
        if (abs(linf - lik[iter - 1]) < tol) {
          conv <- 1
          mess <- paste(G, "group solution converged", sep = "")
          class <- apply(z, 1, function(x) which.max(x))
          s.fin1 <- array(unlist(Sigmain), dim = c(n, n, G))
          p.fin1 <- array(unlist(Psiin), dim = c(p, p, G))
          m.fin <- array(unlist(M), dim = c(n, p, G))
          s.fin2 <- array(0, dim = c(n, n, G))
          p.fin2 <- array(0, dim = c(p, p, G))
          for (e in 1:G) {
            s.fin2[, , e] <- solve(s.fin1[, , e])
            p.fin2[, , e] <- solve(p.fin1[, , e])
          }
          return(list(flag = T, pig = pig, alpha = alpha, eta = eta, M = m.fin, Sigma = s.fin2, Psi = p.fin2, class = class, z = z, v = v, lik = lik[iter]))
        }
      }
    }
    iter <- iter + 1
    if (iter > maxiter) {
      w.message <- paste("Did Not Converge after ", iter - 1, " iterations. Consider increasing max_iter or Tol. G= ", G,
        sep = ""
      )
      return(list(flag = F, exit = w.message))
    }
  }
}
