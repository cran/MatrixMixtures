BIC_Calc <- function(lik, n, p, N, G, mod) {
  if (mod == "MVCN") {
    npar <- G - 1 + G * (n * p + n * (n + 1) / 2 + p * (p + 1) / 2 + 2) - G
  } else if (mod == "MVT") {
    npar <- G - 1 + G * (n * p + n * (n + 1) / 2 + p * (p + 1) / 2 + 1) - G
  } else if (mod == "MVN") {
    npar <- G - 1 + G * (n * p + n * (n + 1) / 2 + p * (p + 1) / 2) - G
  }
  BIC <- 2 * max(lik) - log(N) * npar
}

ConvTestFun <- function(x) {
  if (x$flag == F) {
    x$lik <- c(-Inf)
  }
  return(x)
}
