################################################################################
em.algorithm <- function(Phi, t, v, iters = 100) {
################################################################################
  # Set dimensions
  N <- nrow(Phi)
  M <- ncol(Phi)

  # Set up first iteration
  etas <- (v + 1) / (v + (var(t) ** (-1)) * (t - mean(t)) ** 2 - 2)
  Eta <- diag(etas)  

  # Initialize w and q
  w <- solve(t(Phi) %*% Eta %*% Phi, t(Phi) %*% Eta %*% t)
  q <- ((1 / N) * t(t - Phi %*% w) %*% Eta %*% (t - Phi %*% w)) ** (-1)
  q <- as.numeric(q)

  # All iterations
  for (iter in 1:iters) {
    # Update Eta matrix
    cat('\rIteration:', iter)    
    errors <- t - Phi %*% w  # Error terms
    etas <- (v + 1) / (v + q * (errors ** 2) - 2)
    Eta <- diag(as.numeric(etas))

    # Recompute w and q
    #errors <- t - Phi %*% w  # Error terms
    w <- solve(t(Phi) %*% Eta %*% Phi) %*% t(Phi) %*% Eta %*% t
    q <- ((1 / N) * t(errors) %*% Eta %*% errors) ** (-1)
    q <- as.numeric(q)

    # H matrix
    H <- Eta**(1/2) %*% Phi %*% solve(t(Phi) %*% Eta %*% Phi) %*% t(Phi) %*% Eta**(1/2)
    #H <- Phi %*% solve(t(Phi) %*% Eta %*% Phi) %*% t(Phi) %*% Eta
    #H <- Phi %*% solve(t(Phi) %*% Eta %*% Phi) %*% t(Phi)
    if (iter == iters) { cat('\n') }
  }

  # End
  return(list(w = w, q = q, errors = t - Phi %*% w, etas = etas, H = H))
}

################################################################################
gaussian.model <- function(Phi, t) {
################################################################################
  # Set dimensions
  N <- nrow(Phi)

  # Hat matrix
  H <- Phi %*% solve(t(Phi) %*% Phi) %*% t(Phi)

  # Fitted Values
  t.hat <- H %*% t

  # Observation errors
  errg <- (diag(N) - H) %*% t

  # MLE Estimators
  w.mle <- solve(t(Phi) %*% Phi) %*% t(Phi) %*% t

  # Standard errors for MLE estimators
  q.mle <- ((1 / N) * t(errg) %*% errg) ** (-1)
  var.mle <- solve(as.numeric(q.mle) * t(Phi) %*% Phi)
  sesg <- OpenMx::diag2vec(var.mle) ** (1/2)

  # End
  return(list(w = w.mle, q = q.mle, w.se = sesg, errors = errg, H = H))
}

################################################################################
robust.log.lik <- function(N, t, v, lambda, mu) {
################################################################################
  val <- (N * log(gamma(v / 2 + 1 / 2) / gamma(v / 2)) +
          N * (1 / 2) * log(lambda / (pi * v)) +
          (- (v / 2) - (1 / 2)) * sum(log(1 + (lambda * (t - mu) ** 2) / v)))
  return(as.numeric(val))
}

################################################################################
em.stabilized <- function(Phi, t, v, iters = 100, method) {
################################################################################
  # Set dimensions
  N <- nrow(Phi)
  M <- ncol(Phi)

  # Set up first iteration
  etas <- (v + 1) / (v + (var(t) ** (-1)) * (t - mean(t)) ** 2 - 2)
  Eta <- diag(etas)  

  # Initialize w and q
  w <- solve(t(Phi) %*% Eta %*% Phi, t(Phi) %*% Eta %*% t)
  q <- ((1 / N) * t(t - Phi %*% w) %*% Eta %*% (t - Phi %*% w)) ** (-1)
  q <- as.numeric(q)

  qs <- c()
  ws <- c()
  logliks <- c()
  # All iterations
  for (iter in 1:iters) {
    # Add them up
    ws <- cbind(ws, w)
    qs <- c(qs, q)

    # Update Eta matrix
    cat('\rIteration:', iter)    
    errors <- t - Phi %*% w
    etas <- (v + 1) / (v + q * (errors ** 2) - 2)
    Eta <- diag(as.numeric(etas))

    # Recompute w and q
    w.prev <- w
    q.prev <- q
    errors <- t - Phi %*% w  # Error terms
    w <- solve(t(Phi) %*% Eta %*% Phi) %*% t(Phi) %*% Eta %*% t
    q <- ((1 / N) * t(errors) %*% Eta %*% errors) ** (-1)
    q <- as.numeric(q)
    
    # Log-likelihood  
    mu <- Phi %*% w
    #lambda <- q * ((v - 2) / v)
    lambda <- q * (v / (v - 2))
    #log.lik <- robust.log.lik(N, q, errors, Eta, etas)
    log.lik <- robust.log.lik(N, t, v, lambda, mu)
    logliks <- c(logliks, log.lik)

    # Return
    if (method == 'parameters') {
      if (sum(abs(w - w.prev)) < 1e-6 &&
          sum(abs(q - q.prev)) < 1e-6) {
        cat('\n')
        return(list(ws = ws, qs = qs, n.iter = iter, logliks = logliks))
      }
    } else if (method == 'likelihood') {
      if (iter > 1 && abs(logliks[length(logliks) - 1] - log.lik) < 1e-6) {
        cat('\n')
        return(list(logliks = logliks, n.iter = iter))
      }
    }
    
    # End
    if (iter == iters) {
      cat('\nNo convergence.\n')
    }
  }

  # End
  return(list(ws = ws, qs = qs, n.iter = it))
}


