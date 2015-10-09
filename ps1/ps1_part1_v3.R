################################################################################
# Barcelona Graduate School of Economics
# Master's Degree in Data Science
################################################################################
# Course : Statistical Modelling and Inference
# Title  : Problem Set #1 (part 1)
# Author : (c) Miquel Torrens
# Date   : 2015.10.01
################################################################################
# source('/Users/miquel/Desktop/bgse/courses/term1/smi/ps/ps1/ps1_part1_v3.R')
################################################################################

# Paths
PATH <- '/Users/miquel/Desktop/bgse/courses/term1/smi/'
OUTPUTDIR <- paste(PATH, 'ps/ps1/', sep = '')
DATADIR <- paste(PATH, 'datasets/', sep = '')

# Packages
needed.packages <- c('bbmle', 'OpenMx')
if (length(needed.packages) > 0) {
	for (pkg in needed.packages) {
  	if (! pkg %in% installed.packages()[, 1]) {
  		aux <- try(install.packages(pkg))
  		if (class(aux) != 'try-error') {
	  		cat('Installed package:', pkg, '\n')
  		} else {
  			stop('unable to install package ', pkg)
  		}
		} else {
			aux <- try(library(pkg, character.only = TRUE))
			if (class(aux) != 'try-error') {
				cat('Library:', pkg, '\n')
			} else {
				stop('unable to load package ', pkg)
			}
		}
	}
}

# Load data
file <- paste(DATADIR, 'synthetic_regression.RData', sep = '')
aux <- load(file = file); cat('Loaded file:', file, '\n')
big.data <- get(aux); rm(list = aux[aux != 'big.data']); gc()

# Working data frame
df <- big.data[1:300, 1:31]

################################################################################
# Define estimators
# Phi Matrix
Phi <- as.matrix(cbind(1, df[, 2:ncol(df)]))
colnames(Phi) <- c('X.0', colnames(Phi)[2:ncol(Phi)])

# Hat matrix
H <- Phi %*% solve(t(Phi) %*% Phi) %*% t(Phi)

# Fitted Values
t.hat <- H %*% df[, 1]

# Save results for later exercises
results <- cbind(df[, 1], t.hat)
colnames(results) <- c('t', 't_hat')
save(results, file = paste(OUTPUTDIR, 'model1.RData', sep = '')) 

# Observation errors
err <- (diag(300) - H) %*% df[, 1]

# MLE Estimators
w.mle <- solve(t(Phi) %*% Phi) %*% t(Phi) %*% df[, 1]

# Standard errors for MLE estimators
q.mle <- ((1 / nrow(Phi)) * t(err) %*% err) ** (-1)
var.mle <- solve(as.numeric(q.mle) * t(Phi) %*% Phi)
ses <- OpenMx::diag2vec(var.mle) ** (1/2)

# Relevant parameters: coefficients and standard errors
coefs <- w.mle
se <- ses

################################################################################
# PLOT 1 (WAY 1: plot)
plot(coefs, main = 'MLE Linear Regression Coefficients',
     ylab = 'Coefficient Value', pch = 16, cex = 0.8,
     ylim = c(min(coefs) - 1, max(coefs) + 1),
     xaxt = 'n', xlab = '')

# Print x-axis
axis(side = 1, at = seq_along(coefs), labels = paste('X', 0:30, sep = ''),
     las = 2, cex.axis = 0.7)

# Print upper and lower bound points
points(coefs + qnorm(0.975) * se, pch = 16, cex = 0.5, col = 'red')
points(coefs + qnorm(0.025) * se, pch = 16, cex = 0.5, col = 'red')

# Print lines uniting points
lines(seq_along(coefs), coefs, lty = 'dashed')
lines(seq_along(coefs), coefs + qnorm(0.975) * se, lty = 'dashed', col = 'red')
lines(seq_along(coefs), coefs + qnorm(0.025) * se, lty = 'dashed', col = 'red')

# PLOT 1 (WAY 2: ggplot2)
require(ggplot2)
aux <- data.frame('X_variable' = paste('X', sprintf('%02.0f', 0:30), sep = ''),
                  'Estimated_Coefficient'= coefs,
                  upper = coefs + qnorm(0.975) * se,
                  lower = coefs + qnorm(0.025) * se)
png(paste(OUTPUTDIR, 'ps1_plot1.png', sep = ''), width = 600, height = 400)
ggplot(aux, aes(x = X_variable, y = Estimated_Coefficient)) +
  ggtitle('Confidence intervals for coefficient estimates') +
  theme(plot.title = element_text(lineheight = 0.8, face = "bold")) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymax = upper, ymin = lower))
dev.off()
################################################################################

################################################################################
# PLOT 2
leverage.point <- 3 * length(coefs) / nrow(df)  # Error leverage point
df2 <- data.frame(residuals = ((err - mean(err)) / sd(err)),
                  fitted = t.hat)
#df2[, 3] <- ifelse(df2[, 1] > leverage.point, 2, 1)  # For color
df2[, 3] <- ifelse(diag(H) > leverage.point, 2, 1)  # For color

# Plot
png(paste(OUTPUTDIR, 'ps1_plot2.png', sep = ''), width = 600, height = 400)
plot(df2[, 2], df2[, 1], main = 'Residuals vs. Fitted Values',
     ylab = 'Residual', xlab = 'Fitted value', pch = 16,
     cex = 0.5, ylim = c(min(df2[, 1]) - 0.5,
                         max(df2[, 1]) + 0.5),
     col = df2[, 3])
#abline(h = 3 * length(coefs) / nrow(df), col = 'red', lty = 2)  # Horiz. line
dev.off()
################################################################################

################################################################################
# PLOT 3
set.seed(666)
rand.rsd <- sort(rnorm(300, mean = 0, sd = 1))  # Randomised residuals
#rsd <- sort(err)  # Observed residuals
rsd <- sort(sqrt(q.mle / (1 - diag(H))) * err)

# Plot
png(paste(OUTPUTDIR, 'ps1_plot3.png', sep = ''), width = 600, height = 400)
plot(rsd, rand.rsd, main = 'Standardised vs. Randomised residuals',
     xlab = 'Stadardised residuals',
     ylab = 'Randomised residuals',
     xlim = c(min(rsd), max(rsd)),
     ylim = c(min(rsd), max(rsd)))#, cex = 0.8)

# 45 degree line
abline(0, 1, lty = 2, col = 'red')
dev.off()
################################################################################
