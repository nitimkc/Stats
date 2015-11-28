#source('/Users/miquel/Desktop/bgse/courses/term1/smi/ps/ps6/ps6_sol_v6.R')
# Libraries
library(ggplot2)
library(gridExtra)
library(OpenMx)

# Working directory
setwd('/Users/miquel/Desktop/bgse/courses/term1/smi/ps/ps6')

# Source functions
#source('functionsPS6_v3.R')
source('functionsPS6_v2.R')

# Load the data
aux <- load(file = '../../datasets/synthetic_regression.RData')
data <- get(aux[1])

# Suppress useless stuff
if (! 'data' %in% aux) { rm(list = aux); gc() }

# Choose the data
data <- data[1:300, 1:31]

# Parameters
v <- 10
t <- data[, 1]
Phi <- as.matrix(cbind(1, data[, -1]))
N <- nrow(Phi)
M <- ncol(Phi)
iters <- 100

################################################################################
# EXERCISE 3.1
################################################################################
# Run EM algorithm for Robust model
res <- em.algorithm(Phi, t, v, iters = 100)
w <- res[['w']]
q <- res[['q']]
etas <- res[['etas']]
errors <- res[['errors']]
H <- res[['H']]

# Standard errors for w
A <- matrix(0, M, M)
B <- matrix(0, M, 1)
C <- matrix(0, 1, M)
D <- 0
for (i in 1:N) {
  A <- A + q * (((v + 1) * (v - 2 - q * errors[i] ** 2)) /
               (v + q * errors[i] ** 2 - 2) ** 2) * Phi[i, ] %*% t(Phi[i, ])
  B <- B + q * t(errors[i] ** 3 * ((v + 1) / (v + q * errors[i] ** 2 - 2) ** 2) %*% Phi[i, ])
  C <- C + t(B)
  D <- D + (1 / 2) * sum(1 / (q ** 2) - ((v + 1) * errors[i] ** 4) / (v + q * errors[i] ** 2 - 2) ** 2) 
}

# Build the matrix and obtain the standard errors
Q <- rbind(cbind(A, B), c(C, D))
ses <- sqrt(diag(solve(Q)))[1:M]

# Plot
aux <- data.frame('X_feature' = paste('X', sprintf('%02.0f', 0:30), sep = ''),
                  'Estimated_Coefficient'= w,
                  #upper = w + 1.96 * 10 * ses,
                  #lower = w - 1.96 * 10 * ses)
                  upper = w + 1.96 * ses,
                  lower = w - 1.96 * ses)
plot1 <- ggplot(aux, aes(x = X_feature, y = Estimated_Coefficient)) +
  ggtitle('Robust regression coefficients') +
  theme(plot.title = element_text(lineheight = 0.8, face = 'bold')) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper, ymin = lower))

################################################################################
# Gaussian Model
resG <- gaussian.model(Phi, t)
wG <- resG[['w']]
qG <- as.numeric(resG[['q']])
wG.se <- resG[['w.se']]
errorsG <- resG[['errors']]
HG <- resG[['H']]

# Plot
aux <- data.frame('X_feature' = paste('X', sprintf('%02.0f', 0:30), sep = ''),
                  'Estimated_Coefficient' = wG,
                  upper = wG + 1.96 * wG.se,
                  lower = wG - 1.96 * wG.se)
plot2 <- ggplot(aux, aes(x = X_feature, y = Estimated_Coefficient)) +
  ggtitle('Gaussian Model coefficients') +
  theme(plot.title = element_text(lineheight = 0.8, face = 'bold')) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper, ymin = lower))

png('ps6_plot1_1.png', width = 600, height = 400)
plot1
dev.off()
png('ps6_plot1_2.png', width = 600, height = 400)
plot2
dev.off()

#########################################################################################
# EXERCISE 3.2
#########################################################################################
# Deviance residuals
mu <- Phi %*% w
dr <- rep(NA, length(t))
dg <- rep(NA, length(t))
lambda <- q * (v / (v - 2))
for (n in 1:length(t)) {
  dr[n] <- (-2) * ((log(gamma(v / 2 + 1 / 2) / gamma(v / 2)) +
                   (1 / 2) * log(lambda / (pi * v)) +
                   (- v / 2 - 1 / 2) *
                   log(1 + (lambda * t(t[n]- mu[n]) %*% (t[n] - mu[n])) / v)))
}

# Deviance residuals
devG <- (-1) * log(qG) + qG * errorsG ** 2 + log(2 * pi)
devR <- dr

# Simulating errors
# Robust
set.seed(666)
nsim <- 1000
mu <- rep(0, length(errors))
#Sigma <- q ** (-1) * (diag(length(etas ** (-1))) - H)
Sigma <- q ** (-1) * (diag(length(errors)) - H)
sim <- MASS::mvrnorm(nsim, mu = mu, Sigma = Sigma)
sim <- (-2) * ((log(gamma(v / 2 + 1 / 2) / gamma(v / 2)) +
               (1 / 2) * log(lambda / (pi * v)) +
               (- v / 2 - 1 / 2) *
               log(1 + (lambda * sim ** 2) / v)))

# Gaussian
muG <- rep(0, length(errorsG))
SigmaG <- qG ** (-1) * (diag(length(errorsG)) - HG)
simG <- MASS::mvrnorm(nsim, mu = muG, Sigma = SigmaG)
simG <- (-1) * log(qG) + qG * simG ** 2 + log(2 * pi)

# Cuts
r <- round(quantile(sim, 0.99), 2)
rG <- round(quantile(simG, 0.99), 2)

# Plot
png('ps6_plot2.png', width = 2 * 400, height = 400)
par(mfrow = c(1, 2))
plot(devR, pch = 16, col = 'darkgreen', ylim = c(4, 12),
     main = 'Deviance residuals - Robust regression',
     ylab = 'Deviance residuals')
abline(h = quantile(sim, 0.99), col = 'red')
legend('topright', '99% quantile', col = 'red', lty = 1)
plot(devG, pch = 16, col = 'blue', ylim = c(4, 12),
     main = 'Deviance residuals - Gaussian MLE',
     ylab = 'Deviance residuals')
abline(h = quantile(simG, 0.99), col = 'red')
legend('topright', '99% quantile', col = 'red', lty = 1)
dev.off()

#########################################################################################
# Exercise 3
# OPTION 1: When thetas do not change
stab1 <- em.stabilized(Phi, t, v, iters = 100, method = 'parameters')

# Plot
png('ps6_plot3.png', width = 2 * 600, height = 400)
par(mfrow = c(1, 2))
plot(stab1[['ws']][2, ],
  main = paste('Parameters "w" with EM Algorithm (stabilized at iteration = ',
               stab1[['n.iter']],')', sep = ''),
  xlab = 'Iteration number',
  ylab = 'w_1',
  ylim = c(min(stab1[['ws']][2, ]) - 0.01,
           max(stab1[['ws']][2, ]) + 0.01),
  xaxt = 'n', pch = 16, col = 'darkblue', lty = 2, t = 'b', lwd = 2)
grid()
#abline(h = max(stab2[['logliks']]) + 0.01, col = 'red', lwd = 1)
axis(side = 1, at = seq_along(1:stab1[['n.iter']]),
     labels = 1:stab1[['n.iter']], las = 1, cex.axis = 0.7)

plot(stab1[['qs']],
  main = paste('Parameter "q" with EM Algorithm (stabilized at iteration = ',
               stab1[['n.iter']],')', sep = ''),
  xlab = 'Iteration number',
  ylab = 'q',
  ylim = c(min(stab1[['qs']]) - 0.0005,
           max(stab1[['qs']]) + 0.0005),
  xaxt = 'n', pch = 16, col = 'darkblue', lty = 2, t = 'b', lwd = 2)
grid()
#abline(h = max(stab2[['logliks']]) + 0.01, col = 'red', lwd = 1)
axis(side = 1, at = seq_along(1:stab1[['n.iter']]),
     labels = 1:stab1[['n.iter']], las = 1, cex.axis = 0.7)
dev.off()

# OPTION 2: When MLE does not improve more
stab2 <- em.stabilized(Phi, t, v, iters = 100, method = 'likelihood')

# Plot
png('ps6_plot4.png', width = 600, height = 400)
plot(stab2[['logliks']],
  main = paste('Log-likelihood with the EM Algorithm (stabilized at iteration = ',
               stab2[['n.iter']],')', sep = ''),
  xlab = 'Iteration number',
  ylab = 'Log-likelihood',
  ylim = c(min(stab2[['logliks']]) - 0.2,
           max(stab2[['logliks']]) + 0.2),
  xaxt = 'n', pch = 16, col = 'darkgreen', lty = 2, t = 'b', lwd = 2)
grid()
#abline(h = max(stab2[['logliks']]) + 0.01, col = 'red', lwd = 1)
axis(side = 1, at = seq_along(1:stab2[['n.iter']]),
     labels = 1:stab2[['n.iter']], las = 1, cex.axis = 0.7)
dev.off()

#########################################################################################
# Exercise 4
#nus <- seq(2.5, 100, 0.1)
#nus <- seq(3, 25, 0.1)
nus <- 3:25
#nus <- seq(2.5, 341.5, 1)
opt.nu <- c()
for (nu in nus) {
  out <- em.stabilized(Phi, t, v = nu, iters = 100, method = 'likelihood')
  opt.nu <- rbind(opt.nu, c(nu, out[['logliks']][length(out[['logliks']])]))
}

# When to stop
diffs <- opt.nu[, 2] - c(0, opt.nu[1:(nrow(opt.nu) - 1), 2])
hline <- which(abs(diffs) < 0.1)[1]
stab.nu <- as.numeric(opt.nu[hline, 1])

# Plot
#plot(opt.nu[, 2])
png('ps6_plot5.png', width = 700, height = 400)
plot(opt.nu[, 2],
  main = expression(paste('Choosing ', nu, ' (stabilized at ', nu,
                          #' = ', get(stab.nu), ')', paste = '')),
                          ' = 24)', paste = '')),
  xlab = expression(paste('Value of ', nu)),
  ylab = 'Log-likelihood',
  xaxt = 'n',
  ylim = c(min(opt.nu[, 2]) - 10,
           max(opt.nu[, 2]) + 10),
  pch = 16, col = 'darkgreen', lty = 2, t = 'b', lwd = 2)
grid()
abline(v = hline, col = 'red', lty = 2, lwd = 2)
axis(side = 1, at = 1:23, labels = 3:25, las = 1, cex.axis = 1)
dev.off()
# END OF SCRIPT
