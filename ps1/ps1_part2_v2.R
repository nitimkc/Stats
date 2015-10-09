################################################################################
# Barcelona Graduate School of Economics
# Master's Degree in Data Science
################################################################################
# Course : Statistical Modelling and Inference
# Title  : Problem Set #1 (part 1)
# Author : (c) Miquel Torrens
# Date   : 2015.10.01
################################################################################
# source('/Users/miquel/Desktop/bgse/courses/term1/smi/ps/ps1/ps1_part2_v2.R')
################################################################################

# Paths
PATH <- '/Users/miquel/Desktop/bgse/courses/term1/smi/'
OUTPUTDIR <- paste(PATH, 'ps/ps1/', sep = '')
DATADIR <- paste(PATH, 'datasets/', sep = '')

# Packages
needed.packages <- c('Matrix')
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
df <- big.data[1:300, 1:1001]

# Phi matrix
phi <- as.matrix(cbind(1, df[, 2:ncol(df)]))
colnames(phi) <- c('X.0', colnames(phi)[2:ncol(phi)])

################################################################################
# Singular Value Decomposition
decomp <- svd(t(phi) %*% phi)
#decomp <- svd(phi)

# Singular values
svals <- decomp[['d']]

# Plot them
png(paste(OUTPUTDIR, 'ps1_plot4.png', sep = ''), width = 600, height = 400)
plot(svals[svals > 0], main = 'Singular values of the input matrix',
     ylab = 'Singular values', xlab = 'Index', cex = 0.7)
lines(seq_along(svals[svals > 0]), svals[svals > 0], lty = 'dashed')
dev.off()

png(paste(OUTPUTDIR, 'ps1_plot5.png', sep = ''), width = 600, height = 400)
plot(svals[svals > 1e-10], main = 'Non-zero singular values of the input matrix',
     ylab = 'Singular values', xlab = 'Index', cex = 0.7)
lines(seq_along(svals[svals > 1e-10]), svals[svals > 1e-10], lty = 'dashed')
dev.off()

# Rank of a matrix: number of singular values greater than zero
# https://stat.ethz.ch/pipermail/r-help/2002-February/018865.html
phi.rank <- length(svals[svals > 1e-10])

################################################################################
# Principal Components regression
U <- decomp[['u']][, 1:30]  # U matrix
#U <- decomp[['u']][, 1:31]  # U matrix
Lp <- diag(svals[1:30] ** (-1/2))  # Lambda matrix
#Lp <- diag(svals[1:31] ** (-1/2))  # Lambda matrix
V <- phi %*% U %*% Lp  # V matrix

# Fitted values
t.hat <- V %*% t(V) %*% df[, 1]
#t.hat <- cbind(1, V) %*% rbind(t(V), 1) %*% df[, 1]

########################################################
# Different approach: GO FOR LM WITH SINGULAR VECTORS
pcar <- lm(df[, 1] ~ svd(phi)[['u']][, 1:30])
#pcar <- lm(df[, 1] ~ cbind(1, svd(phi)[['u']][, 1:30]))
t.hat <- fitted(pcar)
########################################################

# Plot fitted vs. observed
txt <- 'Principal Components Regression: Observed versus Fitted values'
png(paste(OUTPUTDIR, 'ps1_plot6.png', sep = ''), width = 600, height = 600)
plot(df[, 1], t.hat,
		 main = txt,
		 ylab = 'Fitted values',
		 xlab = 'Observed values',
		 pch = 16)#, cex = 0.8)
dev.off()

# Previous results
res <- get(load(paste(OUTPUTDIR, 'model1.RData', sep = '')))

# Compare the results graphically
txt1 <- 'MLE Regression: Observed vs. Fitted (30 variables)'
txt2 <- 'PC Regression: Observed vs. Fitted'
png(paste(OUTPUTDIR, 'ps1_plot7.png', sep = ''), width = 1200, height = 600)
par(mfrow = c(1, 2))
plot(res[, 1], res[, 2],
		 main = txt1,
		 ylab = 'Fitted values',
		 xlab = 'Observed values',
		 pch = 16)#, cex = 0.8)
plot(df[, 1], t.hat,
		 main = txt2,
		 ylab = 'Fitted values',
		 xlab = 'Observed values',
		 pch = 16)#, cex = 0.8)
dev.off()

# Compare R^2
r2.orig <- 1 - var(res[, 2] - res[, 1]) / var(res[, 1])
r2.pca <- 1 - var(t.hat - df[, 1]) / var(df[, 1])









