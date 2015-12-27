#source('~/Desktop/bgse/projects/stats/higgs/higgs_v4.R')
# Nedeed packages
require(sgd)
require(ffbase)

# Working directories
PREFIX <- '~/Desktop/bgse/'
DATADIR <- paste(PREFIX, 'datasets/HIGGS/', sep = '')
PATH <- paste(PREFIX, 'projects/stats/higgs/', sep = '')

# Load the Higgs data
file <- paste(PATH, 'higgs500k.RData', sep = '')
if (! file.exists(file)) {
  load.ffdf(paste(DATADIR, 'HIGGSffdf', sep = ''))
  higgs <- Higgs_ffdf[1:500000, ]
  save(higgs, file = file)
} else {
  load(file = file)
}

# Suppress unnecessary features and add define Phi (adding intercept)
higgs[, 'test'] <- NULL
phi <- as.matrix(cbind(1, higgs[, 2:ncol(higgs)]))

# Run normal GLM model
m0 <- glm(signal ~ ., data = higgs, family = 'binomial')
#m0 <- glm(signal ~ . - 1, data = cbind(1, higgs), family = 'binomial')

# Find the weigths
m0w <- m0[['weights']]

# Write TRUE if you want your computer to die
theoretical <- FALSE
if (theoretical) {
  # Uncomputable gamma matrix (500000 x 500000)
  fisher <- t(phi) %*% diag(m0w) %*% phi  # Kill your computer
  vcov <- solve(t(phi) %*% diag(m0w) %*% phi)  # Theoretical approach
  ses <- sqrt(diag(vcov))  # Didn't even get here

  # No help tu use bigalgebra package
  phi <- as.big.matrix(phi)
  gamma <- as.big.matrix(diag(m0w))
  fisher <- t(phi) %*% gamma %*% phi
} else {
  # Instead, we use "optim" function
  y <- higgs[, 'signal']
  x1 <- phi

  # Log-likelihood of the logit
  llk.logit <- function(param, y, x) {
    os <- rep(1,length(x[, 1]))
    x <- cbind(os,x)
    b <- param[1:ncol(x)]
    xb <- x %*% b
    
    # "optim" is a minimizer, so min -ln L(param|y)
    sum(y * log(1 + exp(-xb)) + (1 - y) * log(1 + exp(xb)))
  }

  # Fit logit model using optim
  ls.result <- lm(y ~ x1 - 1)  # Use ls estimates as starting values
  stval <- ls.result$coefficients  # Initial guesses
  
  # Call minimizer procedure (CAREFUL: may take hours!!!)
  file <- paste(PATH, 'logit.RData', sep = '')
  if (! file.exists(file)) {
    logit.m1 <- optim(stval, llk.logit, method = 'BFGS', hessian = TRUE,
                      y = y, x = as.matrix(higgs[, 2:ncol(higgs)]))
    save(logit.m1, file = file)
  } else {
    load(file = file)
  }

  # Get different parameters
  pe.m1 <- logit.m1$par  # Point estimates (SAME AS GLM FUNCTION)
  vc.m1 <- solve(logit.m1$hessian)  # var-cov matrix (what we could not compute)
  hess1 <- logit.m1$hessian  # Hessian matrix
  se.m1 <- sqrt(diag(vc.m1))  # Standard errors (SAME AS GLM FUNCTION)
  ll.m1 <- -logit.m1$value  # Likelihood at maximum

  # Display of the Hessian
  round(hess1, 2)

  # Option 3: simply use "vcov"
  # When I realized about this the rest of the code had already been written
  if (FALSE) {
    mod.IWLS <- glm(signal ~ ., data = higgs, family = binomial())
    beta.mle <- mod.IWLS$coefficients
    beta.cov <- vcov(mod.IWLS)
    beta.sd <- sqrt(diag(beta.cov))
    beta.se <- data.frame(value = beta.mle, se = beta.sd)

    names(beta.mle) <- c('int', names(beta.mle)[2:length(names(beta.mle))])
    names(beta.sd) <- names(beta.mle)
    hess2 <- solve(beta.cov)
    colnames(hess2) <- names(beta.mle)
    rownames(hess2) <- names(beta.mle)
  }
}

################################################################################
# EXERCISE 4
# We generate 50 different orders of the data
set.seed(666)
for (i in 1:50) {
  assign(paste('sam', sprintf('%02.0f', i), sep = ''), sample(1:nrow(higgs)))
}

# Run SGD trials (CAREFUL: this loop takes hours)
file <- paste(DATADIR, 'sgd_results.RData', sep = '')
if (! file.exists(file)) {
  res.df <- c()
  form <- as.formula(signal ~ .)
  for (i in 1:50) {
    s <- get(paste('sam', sprintf('%02.0f', i), sep = ''))
    for (n in c(1, seq(5, 50, 5))) {
      cat('Sample: ', i, ', passes: ', n, '... ', sep = '')
      
      # Run SGD
      res <- sgd(form, data = higgs[s, ], model = 'glm',
                 model.control = list(family = 'binomial'),
                 sgd.control = list(start = pe.m1,
                                    npasses = n,
                                    pass = TRUE,
                                    shuffle = FALSE))
      
      # Store results
      res.df <- cbind(res.df, res$coefficients)
      cat('Done!\n')
    }
  }

  # Save result
  save(res.df, file = file)
} else {
  load(file = file)
}

# Separate coefficients according to number of passes
ses <- c()
for (i in c(1, seq(5, 50, 5))) {
  who <- which(c(1, seq(5, 50, 5)) == i)  
  coefs <- res.df[, seq(who, ncol(res.df), 11)]
  ses <- cbind(ses, apply(coefs, 1, sd))
}

# If we want the first column to be the MLE
ses2 <- ses - se.m1
ses0 <- ses
ses <- cbind(se.m1, ses)

################################################################################
# Plot
png(paste(PATH, 'plot_ex4.png', sep = ''), height = 600, width = 600)
plot(c(0, 0), xlim = c(1, 12), ylim = c(0, 0.03),
     pch = '', xaxt = 'n',
     main = 'Variability of the estimators with SGD',
     xlab = 'Number of passes',
     ylab = 'Standard error of each estimator')
axis(side = 1, at = seq_along(0:11), labels = c('MLE', 1, seq(5, 50, 5)))
grid()

# Add dashed lines for boundaries
abline(v = 1, lty = 2, lwd = 1)
abline(v = 2, lty = 2, lwd = 1)
abline(v = 12, lty = 2, lwd = 1)
abline(h = 0)

# Add legend
legend(x = 'topright', ncol = 2, cex = 0.9, bg = 'white',
       c('MLE', 'SGD') , col = c('red','blue'), lty = 1)

# Add the results for SGD
for (row in 1:nrow(ses)) {
  lines(ses[row, ], col = 'blue')
}

# Add the MLE results
for (se in 1:length(se.m1)) {
  abline(h = se.m1[se], lty = 2, col = 'red', lwd = 1)
}
dev.off()

png(paste(PATH, 'plot_ex4_v2.png', sep = ''), height = 600, width = 600)
plot(c(0, 0), xlim = c(1, 11), ylim = c(-0.019, 0.006),
     pch = '', xaxt = 'n',
     main = 'Variablity of estimators: MLE minus SGD',
     xlab = 'Number of passes',
     ylab = 'Gap of the standard error of each estimator')
axis(side = 1, at = seq_along(1:11), labels = c(1, seq(5, 50, 5)))
grid()

# Add dashed lines for boundaries
abline(v = 1, lty = 2, lwd = 1)
abline(v = 11, lty = 2, lwd = 1)
abline(h = 0)

# Add the results for SGD
for (row in 1:nrow(ses)) {
  lines(ses2[row, ], col = 'darkgreen')
}
dev.off()

png(paste(PATH, 'plot_ex4_v3.png', sep = ''), height = 700, width = 700)
plot(c(1, 0), xlim = c(1, 12), ylim = c(0, 0.03),
     pch = '', xaxt = 'n',
     main = 'Variability of the estimators with SGD',
     xlab = 'Number of passes',
     ylab = 'Standard error of each estimator')
axis(side = 1, at = seq_along(0:11), labels = c('MLE', 1, seq(5, 50, 5)))
grid()

# Add dashed lines for boundaries
abline(v = 1, lty = 2, lwd = 1)
abline(v = 2, lty = 2, lwd = 1)
abline(v = 12, lty = 2, lwd = 1)
abline(h = 0)

# Add legend
legend(x = 'topright', ncol = 2, cex = 0.9, bg = 'white',
       c('MLE', 'SGD') , col = c('red','blue'), lty = 1)

# Add the results for SGD
for (row in 1:nrow(ses)) {
  lines(ses[row, ], col = 'blue')
  lines(ses[row, 1:2], col = 'red', lwd = 2)
  lines(ses[row, 1:2], col = 'white', lty = 2, lwd = 2)
}

# Add MLE points
points(cbind(1, ses[, 1]), pch = 16, col = 'red')
dev.off()
# END OF SCRIPT
