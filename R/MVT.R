pig_updateT <- function(z, N) {
  pig <- colSums(z) / N
  return(pig)
}

Mg_UpdateT <- function(X, z, estep) {
  Mg <- list()
  G <- ncol(z)
  N <- length(X)
  stemp <- z * estep[[1]]
  for (g in 1:G) {
    numer <- Reduce("+", lapply(1:N, function(i) stemp[i, g] * X[[i]]))
    Mg[[g]] <- numer / sum(stemp[, g])
  }
  return(Mg)
}

Sigma_UpdateT <- function(X, z, w, M, Psiin) {
  Sigmain <- list()
  p <- dim(Psiin[[1]])[1]
  detSig <- NULL
  G <- ncol(z)
  N <- nrow(z)
  stemp <- z * w
  for (g in 1:G) {
    numer <- Reduce("+", lapply(1:N, function(i) stemp[i, g] * (X[[i]] - M[[g]]) %*% Psiin[[g]] %*% t(X[[i]] - M[[g]])))
    temp <- numer / (p * colSums(z)[g])
    Sigmain[[g]] <- solve(temp)
    detSig[g] <- det(Sigmain[[g]])
  }
  return(list(Sigmain, detSig, prob = 0))
}

Psi_UpdateT <- function(X, z, w, M, Sigmain) {
  Psiin <- list()
  n <- dim(Sigmain[[1]])[1]
  detPsi <- NULL
  G <- ncol(z)
  N <- nrow(z)
  for (g in 1:G) {
    stemp <- z * w
    numer <- Reduce("+", lapply(1:N, function(i) stemp[i, g] * t(X[[i]] - M[[g]]) %*% Sigmain[[g]] %*% (X[[i]] - M[[g]])))
    temp <- numer / (n * colSums(z)[g])
    Psiin[[g]] <- solve(temp)
    detPsi[g] <- det(Psiin[[g]])
  }
  return(list(Psiin, detPsi, prob = 0))
}

nug_updateT <- function(estep, z) {
  nu <- NULL
  G <- ncol(z)
  for (g in 1:G) {
    Ng <- sum(z[, g])
    dfnew <- try(stats::uniroot(function(v) log(v / 2) + 1 - digamma(v / 2) + (1 / Ng) * sum(z[, g] * (estep[[2]][, g] - estep[[1]][, g])), lower = 0.001, upper = 1000)$root, silent = F)
    if (is.double(dfnew)) {
      if (dfnew > 200) {
        nu[g] <- 200
      } else if (dfnew < 2) {
        nu[g] <- 2
      } else {
        nu[g] <- dfnew
      }
    } else {
      w.message <- "Difficulty Updating Degrees of Freedom"
      return(list(w.message, prob = 1))
    }
  }
  return(list(nu = nu, prob = 0))
}

MVt <- function(X, M, Sigmain, Psiin, nu, n, p, N) {
  lambda <- -(nu + n * p) / 2
  Delta <- sum(diag(Sigmain %*% (X - M) %*% Psiin %*% t(X - M))) + nu
  dens <- log(det(Sigmain)) * (p / 2) + log(det(Psiin)) * (n / 2) + log(gamma((n * p + nu) / 2)) - (n * p / 2) * (log(pi) + log(nu)) - log(gamma(nu / 2)) - ((n * p + nu) / 2) * (log(Delta) - log(nu))
  return(list(dens = exp(dens)))
}

zig_updateT <- function(X, M, Sigmain, Psiin, N, nu, n, p, pig, G) {
  zmat <- matrix(NA, nrow = N, ncol = G)
  for (g in 1:G) {
    zmat[, g] <- pig[g] * unlist(lapply(X, function(x) MVt(x, M[[g]], Sigmain[[g]], Psiin[[g]], nu[g], n, p, N)$dens))
  }
  if (any(is.na(zmat))) {
    w.message <- "NA in Zig update!"
    return(list(w.message, prob = 1))
  }
  return(list(zmat = zmat, prob = 0))
}

Estep2T <- function(X, M, Sigmain, Psiin, nu) {
  G <- length(M)
  N <- length(X)
  n <- dim(Sigmain[[1]])[1]
  p <- dim(Psiin[[1]])[1]
  wupdate <- matrix(0, N, G)
  logwupdate <- matrix(0, N, G)
  for (g in 1:G) {
    Delta <- unlist(lapply(X, function(x) sum(diag(Sigmain[[g]] %*% (x - M[[g]]) %*% Psiin[[g]] %*% t(x - M[[g]]))) + nu[g]))
    wupdate[, g] <- (n * p + nu[g]) / Delta
    logwupdate[, g] <- digamma(((n * p + nu[g]) / 2)) - log(0.5 * Delta)
  }
  return(list(wupdate, logwupdate))
}

EM_MVT <- function(X, G, maxiter = 1000, tol = 1.0e-03, sm = 1, init = "rand") {
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
  for (g in 1:G) {
    M[[g]] <- Reduce("+", lapply(1:N, function(i) z[i, g] * X[[i]])) / sum(z[, g])
    Sigma <- Reduce("+", lapply(1:N, function(i) z[i, g] * (X[[i]] - M[[g]]) %*% t(X[[i]] - M[[g]]))) / (p * sum(z[, g]))
    Sigmain[[g]] <- solve(Sigma)
    detSig[g] <- det(Sigmain[[g]])
    Psi <- Reduce("+", lapply(1:N, function(i) z[i, g] * t(X[[i]] - M[[g]]) %*% (X[[i]] - M[[g]]))) / (n * sum(z[, g]))
    Psiin[[g]] <- solve(Psi)
    detPsi[g] <- det(Psiin[[g]])
    nu <- rep(50, G)
  }
  estep <- Estep2T(X, M, Sigmain, Psiin, nu)
  conv <- 0
  prob <- 0
  iter <- 1
  lik <- NULL
  while (conv == 0) {
    pig <- pig_updateT(z, N)
    M <- Mg_UpdateT(X, z, estep)
    nutemp <- nug_updateT(estep, z)
    if (nutemp$prob == 1) {
      return(list(flag = 1))
    }
    nu <- nutemp[[1]]
    Sigmatest <- Sigma_UpdateT(X, z, estep[[1]], M, Psiin)
    if (Sigmatest$prob == 1) {
      return(list(flag = F))
    }
    Sigmain <- Sigmatest[[1]]
    detSig <- Sigmatest[[2]]
    Psitest <- Psi_UpdateT(X, z, estep[[1]], M, Sigmain)
    if (Psitest$prob == 1) {
      return(list(flag = F))
    }
    Psiin <- Psitest[[1]]
    detPsi <- Psitest[[2]]
    Estep <- zig_updateT(X, M, Sigmain, Psiin, N, nu, n, p, pig, G)
    zmat <- Estep[[1]]
    z <- zmat / rowSums(zmat)
    estep <- Estep2T(X, M, Sigmain, Psiin, nu)
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
          return(list(flag = T, pig = pig, nu = nu, M = m.fin, Sigma = s.fin2, Psi = p.fin2, class = class, z = z, lik = lik[iter]))
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
