################################################################################
new.trial <- function(g, delta = 0.5) {
################################################################################
  #delta <- 0.5  # This is not explicitely defined
  a <- (-1) * (1 / 2)
  b <- 0
  D <- delta * diag(nrow(phi))

  # K matrix
  K <- g * phi %*% solve(D + g * t(phi) %*% phi) %*% t(phi)

  # New parameters
  a. <- a + nrow(phi) / 2
  b. <- b + (1 / 2) *  t(t) %*% (diag(nrow(K)) - K) %*% t
  w.bayes <- g * solve(D + g * t(phi) %*% phi) %*% t(phi) %*% t
  D. <- g ** (-1) * (D + g * t(phi) %*% phi)

  # Return
  return(var(t - (phi %*% w.bayes)) ** (-1))
}

################################################################################
# Find optimal "g"
find.g <- function(delta = 1) {
################################################################################
  g0 <- var(t - mean(t)) ** (-1)  # Model with just the intercept
  res <- g0
  res[2] <- new.trial(g = g0, delta)
  for (i in 3:100) {
    cat('\rTrial:', i)
    res[i] <- new.trial(g = res[i - 1], delta)
    if (i > 5 && round(res[i], 5) == round(res[i - 5], 5)) {
      cat('\ng has stabilized at:', res[i], '\n')
      return(invisible(res[i]))
    }
  }
}

################################################################################
plot.model <- function(delta, dpng = NULL) {
################################################################################
  # Parameters
  #delta <- 0.5  # This is not explicitely defined
  g <- find.g(delta = delta)
  #g <- 10       # This is not explicitely defined
  a <- (-1) * (1 / 2)
  b <- 0
  D <- delta * diag(nrow(phi))

  # K matrix
  K <- g * phi %*% solve(D + g * t(phi) %*% phi) %*% t(phi)

  # New parameters
  a. <- a + nrow(phi) / 2
  b. <- b + (1 / 2) *  t(t) %*% (diag(nrow(K)) - K) %*% t
  w.bayes <- g * solve(D + g * t(phi) %*% phi) %*% t(phi) %*% t
  D. <- g ** (-1) * (D + g * t(phi) %*% phi)

  # Replicate Phi with more observations
  M <- 9
  new.x <- seq(0, 1, 1 / 10000)
  phi3 <- matrix(nrow = length(new.x), ncol = M + 1)
  invisible(sapply(1:length(new.x), function(i) {
    phi3[i, ] <<- phix(x = new.x[i], M = 9, basis = 'Gauss')
  }))

  # Add new phi?
  par1 <- phi3 %*% w.bayes
  par2 <- rep(0, length = nrow(phi3))
  invisible(sapply(1:nrow(phi3), function(i) {
    par2[i] <<- (a. / b.) * 1 / (1 + t(phi3[i, ]) %*% D. %*% phi3[i, ])
    #par2[i] <<- 1 / (1 + t(phi3[i, ]) %*% D. %*% phi3[i, ])
  }))
  par3 <- 2 * a.

  par2. <- par2
  #par2. <- sqrt(((par3 / (par3 - 2)) * (1 / ((a. / b.) * par2))))
  #par2. <- sqrt(diag((a. / (a. - 2)) * solve(D.)))
  #par2. <- 1 / par2

  if (! is.null(dpng)) {
    png(paste0('~/Desktop/bgse/courses/term1/smi/ps/ps3/', dpng))
  }
  main <- paste('Bayesian linear predictor (delta = ', delta, ', approx. g = ',
                round(g, 0), ')', sep = '')
  plot(x, t, ylim = c(min(t) - 1.5, max(t) + 1.5), pch = 16,
       main = main, xlab = 'Input', ylab = 'Output')
  lines(new.x, par1, col = 'blue')
  lines(new.x, par1 + sqrt(par2.), col = 'red')
  lines(new.x, par1 - sqrt(par2.), col = 'red')
  legend('topright', c('predicted mean', '+/- standard deviation'),
         col = c('blue', 'red'), lty = 1)
  if (! is.null(dpng)) {
    dev.off()
  }
    
  return(invisible(list(par1 = par1, par2 = par2., par3 = par3,
                        w.bayes = w.bayes, D. = D., a. = a., b. = b.)))
}