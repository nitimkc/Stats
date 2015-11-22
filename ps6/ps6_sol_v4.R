
# Libraries
library(ggplot2)
library(gridExtra)
library(OpenMx)

# Working directory
setwd('/Users/miquel/Desktop/bgse/courses/term1/smi/ps/ps6')

# Source functions
source('functionsPS6.R')

# Load the data
aux <- load(file = '../../datasets/synthetic_regression.RData')
data <- get(aux[1])

# Suppress useless stuff
if (! 'data' %in% aux) { rm(list = aux); gc() }

# Choose the data
data <- data[1:300, 1:31]

# Save data
save(data, file = '/Users/miquel/Desktop/data.RData')

# Parameters
v <- 10
Phi <- cbind(1, data)
iters <- 100
t <- Phi[, 2]
Phi <- as.matrix(Phi[, -2])

################################################################################
# EXERCISE 3.1
################################################################################
# Run EM algorithm for Robust model
res <- em.algorithm(Phi, t, v, iters = 100)
q <- res[['q']]
etas <- res[['etas']]
errors <- res[['errors']]

# Standard errors for w
el1 <- q * sum(((v + 1) * (v - 2 - q * errors ** 2)) /
                (v + q * errors ** 2 - 2) ** 2) * t(Phi) %*% Phi
el2 <- q * t(err^3 * ((v + 1) / (v + q * errors ** 2 - 2) ** 2)) %*% Phi
el3 <- t(el2)
el4 <- (1 / 2) * sum((1 / (q ** 2)) - ((v + 1) * errors ** 4) /
                                       (v + q * errors ** 2 - 2) ** 2)

# Build the matrix and obtain the standard errors
Q <- rbind(cbind(el1, t(el2)), c(t(el3), el4))
#ses <- sqrt(diag(solve(Q)))[1:31]
ses <- sqrt(diag(solve(Q)))[1:M]

# Plot
aux <- data.frame('X_feature' = paste('X', sprintf('%02.0f', 0:30), sep = ''),
                  'Estimated_Coefficient'= w,
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

png('ps6_plot1.png', width = 2 * 600, height = 400)
grid.arrange(plot1, plot2, ncol = 2)
dev.off()

#########################################################################################
# EXERCISE 3.2
#########################################################################################
# Deviance residuals
devG <- - N * log(qG) + qG * errorsG ** 2 + 2 * log(2 * pi)
devR <- - N * log(q) + q * errors ** 2 + 2 * log(2 * pi) - N * log(diag(Eta)) + diag(Eta)

par(mfrow = c(1, 2))
plot(devR / sum(devR), pch = 16, col = 'darkgreen', ylim = c(0.0032, 0.0045))
abline(h = quantile(devR / sum(devR), 0.99), col = 'red')
plot(devG / sum(devG), pch = 16, col = 'blue', ylim = c(0.0032, 0.0045))
abline(h = quantile(devG / sum(devG), 0.99), col = 'red')

#devG <- - N * log(qG) + qG * errorsG ** 2
#devR <- - N * log(q) + q * errors ** 2 - N * log(diag(Eta)) + diag(Eta)

pred <- Phi %*% w
sign(err) * (2 * pred * log(t / pred) + 2 * (N - t) * log((N - t) / (N - pred))) ** (1/2)


dev.gauss <- q.mle * (errg ** 2)
devR <- q[1, 1] * (err ** 2)# %*% diag(Eta)


devG <- errorsG ** 2
devR <- (errors ** 2) * etas



png('ps6_plot2.png', width = 2 * 400, height = 400)
#png('ps6_plot2.png', width = 2 * 600, height = 400)
par(mfrow = c(1, 2))
plot(devR, pch = 16, col = 'darkgreen', ylim = c(0, max(devG) + 5),
     main = 'Deviance residuals - Robust regression',
     ylab = 'Deviance residuals')
#plot(devR, pch = 16, col = 'darkgreen')
#abline(h = quantile(devR, 0.99), col = 'red')
abline(h = quantile(devR, 0.99), col = 'red')
legend('topright', '99% quantile', col = 'red', lty = 1)
plot(devG, pch = 16, col = 'blue', ylim = c(0, max(devG) + 5),
     main = 'Deviance residuals - Gaussian MLE',
     ylab = 'Deviance residuals')
#plot(devG, pch = 16, col = 'blue')
abline(h = quantile(devG, 0.99), col = 'red')
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
  main = paste('Parameters "w" with the EM Algorithm (stabilized at iteration = ',
               stab1[['n.iter']],')', sep = ''),
  xlab = 'Iteration number',
  ylab = 'phi_1',
  ylim = c(min(stab1[['ws']][2, ]) - 0.01,
           max(stab1[['ws']][2, ]) + 0.01),
  xaxt = 'n', pch = 16, col = 'darkblue', lty = 2, t = 'b', lwd = 2)
grid()
#abline(h = max(stab2[['logliks']]) + 0.01, col = 'red', lwd = 1)
axis(side = 1, at = seq_along(1:stab1[['n.iter']]),
     labels = 1:stab1[['n.iter']], las = 1, cex.axis = 0.7)

plot(stab1[['qs']],
  main = paste('Parameter "q" with the EM Algorithm (stabilized at iteration = ',
               stab1[['n.iter']],')', sep = ''),
  xlab = 'Iteration number',
  ylab = 'q',
  ylim = c(min(stab1[['qs']]) - 0.001,
           max(stab1[['qs']]) + 0.001),
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
nus <- seq(2.5, 100, 0.1)
#nus <- seq(0.1, 100, 0.1)
opt.nu <- c()
for (nu in nus) {
  out <- em.stabilized(Phi, t, v = nu, iters = 100, method = 'likelihood')
  opt.nu <- rbind(opt.nu, c(nu, max(out[[1]])))
}

# Plot
#plot(opt.nu[, 2])
png('ps6_plot5.png', width = 600, height = 400)
plot(opt.nu[, 2],
  main = 'Choosing "nu"',         
  xlab = 'Value of "nu"',
  ylab = 'Log-likelihood',
  ylim = c(min(opt.nu[, 2]) - 0.2,
           max(opt.nu[, 2]) + 0.2),
  pch = 16, col = 'darkgreen', lty = 2, t = 'b', lwd = 2)
  #xaxt = 'n', pch = 16, col = 'darkgreen', lty = 2, t = 'b', lwd = 2)
grid()
#axis(side = 1, at = seq_along(opt.nu[, 1]),
#     labels = opt.nu[, 1], las = 1, cex.axis = 0.7)
dev.off()










