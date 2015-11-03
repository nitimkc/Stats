################################################################################
# Barcelona Graduate School of Economics
# Master's Degree in Data Science
################################################################################
# Course : Statistical Modelling and Inference
# Title  : Problem Set #4
# Author : (c) Miquel Torrens
# Date   : 2015.11.01
################################################################################
# source('/Users/miquel/Desktop/bgse/courses/term1/smi/ps/ps4/smi_ps4_v1.R')
################################################################################

# Path
PATH <- '/Users/miquel/Desktop/bgse/courses/term1/smi/ps/ps4/'
DATADIR <- paste(PATH, 'ARM_Data/', sep = '')

# Widths and names
widths <- c(5, 1, 2, 3, 3, 1, 1, 2, 2, rep(1, 20), rep(2, 12), 1, 3, 1,
	        rep(3, 10), 1, 2, 1, 1, 1, 1, 1, 2, 3, 3, 1, 3, rep(1, 19), 2, 1, 1,
	        1, 3, 3, rep(1, 40), 2, 2, 2, 2, 1, 1, 1, 1, 1, 6, 2, 6, 2, 1)
cols <- c('ID', 'EMP', 'EMPOTH', 'OCC', 'IND', 'SELFEMP', 'ROUTINE', 'TASK1',
	      'TASK2', 'ROUTINE2', 'GOOD', 'OTHGOOD', 'THANK', 'WENJOY', 'LEARN',
	      'RECOG', 'WORKSAT', 'RECOM', 'JOBHOME', 'SUP', 'DECHOW', 'DECWHAT',
	      'DISAG', 'PROMOTE', 'SUPERV', 'SUPERV2', 'GOALS', 'MANAG', 'MANLEV',
	      'ADULTS', 'KIDS', 'AGEKID1', 'AGEKID2', 'AGEKID3', 'AGEKID4',
	      'AGEKID5', 'AGEKID6', 'AGEKID7', 'AGEKID8', 'KIDCARE', 'KIDCOTH',
	      'DIFCARE', 'MONCARE', 'STRNCARE', 'COOK', 'SHOP', 'CLEAN', 'LAUNDRY',
	      'REPAIR', 'DISHES', 'BUDGET', 'PLANS', 'CHILDC', 'HSWORK', 'MARSTAT',
	      'YRALONE', 'PARTNER', 'MARHAPPY', 'CHANGE', 'DIVTHOT', 'SPEMP',
	      'SPEMPOTH', 'SPOCC', 'SPIND', 'SPFEEL', 'SPHSWORK', 'VAC', 'HOUSE',
	      'MOVE', 'BUY', 'STRNMED', 'STRNFOOD', 'STRNBILL', 'WORRY', 'TENSE',
	      'RESTLESS', 'AFRAID', 'FEAR', 'MAD', 'YELL', 'ANGRY', 'TRUST', 'SUSP',
	      'AGAINST', 'HEALTH', 'WALK', 'FARWALK', 'EXER', 'DIET', 'HEIGHT',
	      'WEIGHT', 'SMOKENOW', 'SMOKEV', 'STAIRS', 'KNEEL', 'CARRY', 'HAND',
	      'SEE', 'HEAR', 'DIFWALK', 'PAIN', 'HEAD', 'WEAK', 'SLEEP', 'EFFORT',
	      'GETGO', 'MIND', 'SAD', 'LONELY', 'BLUE', 'ENJOY', 'HOPE', 'HAPPY',
	      'FATGOOD', 'FATHAPPY', 'RESPSUC', 'RESPANY', 'FATPROB', 'FATBAD',
	      'RESPMIS', 'RESPFAIL', 'EMOT', 'SUPTURN', 'SUPTALK', 'USGOODL',
	      'USACHIEV', 'USDES', 'USEFFORT', 'USBADL', 'USGREED', 'OWN', 'ED',
	      'MOMED', 'FATHED', 'YEARBN', 'RACE', 'RACEOTH', 'HISP', 'REL',
	      'RELOTH', 'EARN1', 'EARN2', 'FAMINC1', 'FAMINC2', 'SEX')

# Load data
wfw <- read.fwf(paste(PATH, 'wfw90.dat', sep = ''), widths)

# Set names
colnames(wfw) <- cols

#library('foreign')
#heights <- read.dta(paste(PATH, 'heights.dta', sep = ''))

################################################################################
# EXERCISE 2
# EXERCISE 2.1 #################################################################
dta <- wfw[, c('EARN1', 'EARN2', 'SEX', 'HEIGHT', 'WEIGHT')]

# Define total earnings
dta[, 'EARNT'] <- apply(dta, 1, function(x) {
  #sum(x['EARN1'], 1e3 * x['EARN2'], na.rm = TRUE)
  sum(x['EARN1'] / 1e3, x['EARN2'], na.rm = TRUE)
})

# Indicator on whether the earnings are exact or not
dta[, 'INEXACT'] <- 1
dta[which(! is.na(dta[, 'EARN1'])), 'INEXACT'] <- 0

# Set height in inches
dta[, 'HEIGHT_I'] <- 12 * as.numeric(substr(dta[, 'HEIGHT'], 1, 1)) +
                     as.numeric(substr(dta[, 'HEIGHT'], 2, 3)) 

# Set sex as 0/1 binary
dta[, 'MEN'] <- ifelse(dta[, 'SEX'] == 1, 1, 0)

# Clean data set of extreme values
dta <- dta[which(dta[, 'WEIGHT'] < 500), ]
dta <- dta[which(dta[, 'HEIGHT'] < 800), ]
dta <- dta[which(dta[, 'EARNT'] < 300), ]

# EXERCISE 2.2 #################################################################
# Define centered variables
dta[, 'EARNT_C'] <- dta[, 'EARNT'] - mean(dta[, 'EARNT'])
dta[, 'WEIGHT_C'] <- dta[, 'WEIGHT'] - mean(dta[, 'WEIGHT'])
dta[, 'HEIGHT_I_C'] <- dta[, 'HEIGHT_I'] - mean(dta[, 'HEIGHT_I'])

# Standardised variables
dta[, 'HEIGHT_I_ST'] <- (dta[, 'HEIGHT_I'] - mean(dta[, 'HEIGHT_I'])) / sd(dta[, 'HEIGHT_I'])
dta[, 'WEIGHT_ST'] <- (dta[, 'WEIGHT'] - mean(dta[, 'WEIGHT'])) / sd(dta[, 'WEIGHT'])

# Linear model
m00 <- lm(EARNT ~ HEIGHT_I + INEXACT, data = dta)
m01 <- lm(EARNT ~ HEIGHT_I, data = dta)
m02 <- lm(EARNT ~ HEIGHT_I_C, data = dta)
#m02 <- lm(EARNT_C ~ HEIGHT_I_C + EXACT, data = dta)

# EXERCISE 2.3 #################################################################
dta2 <- dta[which(dta[, 'EARNT'] > 0), ]

# Linear-linear
m03 <- lm(EARNT ~ HEIGHT_I + WEIGHT + MEN + INEXACT, data = dta)
m04 <- lm(EARNT ~ HEIGHT_I + WEIGHT + MEN, data = dta)
m05 <- lm(EARNT ~ HEIGHT_I + WEIGHT + MEN + WEIGHT * MEN + HEIGHT_I * MEN + INEXACT, data = dta)
m06 <- lm(EARNT ~ HEIGHT_I + WEIGHT + MEN + WEIGHT * MEN + HEIGHT_I * MEN, data = dta)

# Log-linear
m07 <- lm(log(EARNT) ~ HEIGHT_I + WEIGHT + MEN + WEIGHT * MEN, data = dta2)
m08 <- lm(log(EARNT) ~ HEIGHT_I + WEIGHT + MEN + WEIGHT * MEN + INEXACT, data = dta2)
m09 <- lm(log(EARNT) ~ HEIGHT_I + WEIGHT + MEN + WEIGHT * MEN + HEIGHT * MEN, data = dta2)
m10 <- lm(log(EARNT) ~ HEIGHT_I + WEIGHT + MEN + WEIGHT * MEN, data = dta2)
m11 <- lm(log(EARNT) ~ HEIGHT_I_ST + WEIGHT_ST + MEN + WEIGHT_ST * MEN + INEXACT, data = dta2)
m12 <- lm(log(EARNT) ~ HEIGHT_I_ST + WEIGHT_ST + MEN + WEIGHT_ST * MEN + HEIGHT_I_ST * MEN + INEXACT, data = dta2)
m13 <- lm(log(EARNT) ~ HEIGHT_I_ST + WEIGHT_ST + MEN + WEIGHT_ST * MEN + HEIGHT_I_ST * MEN + INEXACT * MEN, data = dta2) # !!!
#m13 <- lm(log(EARNT) ~ HEIGHT_I_ST + WEIGHT_ST + MEN + WEIGHT_ST * MEN + HEIGHT_I_ST * MEN + INEXACT, data = dta2)
m14 <- lm(log(EARNT) ~ HEIGHT_I_ST + WEIGHT_ST + MEN + WEIGHT_ST * MEN + HEIGHT_I_ST * MEN, data = dta2)
m15 <- lm(log(EARNT) ~ HEIGHT_I + WEIGHT + MEN + WEIGHT * MEN + HEIGHT_I * MEN + INEXACT, data = dta2)
m16 <- lm(log(EARNT) ~ HEIGHT_I + WEIGHT + MEN + WEIGHT * MEN + INEXACT, data = dta2)
m17 <- lm(log(EARNT) ~ HEIGHT_I + WEIGHT + MEN + INEXACT, data = dta2)
m18 <- lm(log(EARNT) ~ HEIGHT_I + WEIGHT + MEN, data = dta2)
m19 <- lm(log(EARNT) ~ HEIGHT_I + WEIGHT + MEN + WEIGHT * MEN + HEIGHT_I * MEN + INEXACT, data = dta2)
m20 <- lm(log(EARNT) ~ HEIGHT_I + WEIGHT + MEN + WEIGHT * MEN + HEIGHT_I * MEN, data = dta2)

# Log-log
m21 <- lm(log(EARNT) ~ log(HEIGHT_I) + log(WEIGHT) + MEN + INEXACT, data = dta2)

# Evaluate R-squared
for (i in sprintf('%02.0f', 1:21)) {
  mod <- get(paste0('m', i))
  cat('m', i, ': ', round(summary(mod)$adj.r.squared, 4), ' resid: ',
  		round(sqrt(deviance(mod) / df.residual(mod)), 2), '\n', sep = '')
}  
################################################################################

################################################################################
# EXERCISE 3
# EXERCISE 3.1 #################################################################
set.seed(666)
p1 <- rnorm(n = 1000, mean = 600, sd = 400)
p2 <- rnorm(n = 1000, mean = 3L, sd = 1L)

# Plot the simulated draws
png(paste(PATH, 'plot_ex3dot1.png', sep = ''))
plot(p1, p2, pch = 16,# col = 'gray',
	 main = 'Simulation of 1,000 draws',
	 xlab = 'Cost difference',
	 ylab = 'Effectiveness difference')
abline(v = 0, lty = 'dashed', col = 'red', lwd = 2)
abline(h = 0, lty = 'dashed', col = 'red', lwd = 2)
dev.off()

# EXERCISE 3.2 #################################################################
# 50% interval
qnorm(0.25, mean = 600, sd = 400) / qnorm(0.25, mean = 3L, sd = 1L)
qnorm(0.75, mean = 600, sd = 400) / qnorm(0.75, mean = 3L, sd = 1L)

# 95% interval
qnorm(0.025, mean = 600, sd = 400) / qnorm(0.025, mean = 3L, sd = 1L)
qnorm(0.975, mean = 600, sd = 400) / qnorm(0.975, mean = 3L, sd = 1L)

# Altogether
if (FALSE) {
	qnorm(0.25, mean = 600 / 3L, sd = 400 / 1L)
	qnorm(0.75, mean = 600 / 3L, sd = 400 / 1L)
	qnorm(0.025, mean = 600 / 3L, sd = 400 / 1L)
	qnorm(0.975, mean = 600 / 3L, sd = 400 / 1L)
}

# Using quantiles
q.ratio1 <- quantile(p1 / p2, c(0.025, 0.25, 0.75, 0.975))
#t.test(p1 / p2)

# EXERCISE 3.3 #################################################################
set.seed(666)
p1 <- rnorm(n = 1000, mean = 600, sd = 400)
p2 <- rnorm(n = 1000, mean = 3L, sd = 2L)

# Plot the simulated draws
png(paste(PATH, 'plot_ex3dot3.png', sep = ''))
plot(p1, p2, pch = 16,# col = 'gray',
	 main = 'Simulation of 1,000 draws',
	 xlab = 'Cost difference',
	 ylab = 'Effectiveness difference')
abline(v = 0, lty = 'dashed', col = 'red', lwd = 2)
abline(h = 0, lty = 'dashed', col = 'red', lwd = 2)
dev.off()

# 50% interval
qnorm(0.25, mean = 600, sd = 400) / qnorm(0.25, mean = 3L, sd = 2L)
qnorm(0.75, mean = 600, sd = 400) / qnorm(0.75, mean = 3L, sd = 2L)

# 95% interval
qnorm(0.025, mean = 600, sd = 400) / qnorm(0.025, mean = 3L, sd = 2L)
qnorm(0.975, mean = 600, sd = 400) / qnorm(0.975, mean = 3L, sd = 2L)

# Using quantiles
q.ratio2 <- quantile(p1 / p2, c(0.025, 0.25, 0.75, 0.975))

################################################################################
# EXERCISE 4
# Read data
congress <- vector('list', 49)
for (i in 1:49) {
  year <- 1896 + 2 * (i - 1)
  file <- paste(DATADIR, 'cong3/', year, '.asc', sep = '')
  data.year <- matrix(scan(file), byrow = TRUE, ncol = 5)
  data.year <- cbind(rep(year, nrow(data.year)), data.year)
  congress[[i]] <- data.year
}

# Tidy data from period 1986-1990
i86 <- (1986 - 1896) / 2 + 1
cong86 <- congress[[i86]]
cong88 <- congress[[i86 + 1]]
cong90 <- congress[[i86 + 2]]

v86 <- cong86[, 5] / (cong86[, 5] + cong86[, 6])
bad86 <- cong86[, 5] == -9 | cong86[, 6] == -9
v86[bad86] <- NA
contested86 <- v86 > 0.1 & v86 < 0.9
inc86 <- cong86[, 4]

v88 <- cong88[, 5] / (cong88[, 5] + cong88[, 6])
bad88 <- cong88[, 5] == -9 | cong88[, 6] == -9
v88[bad88] <- NA
contested88 <- v88 > 0.1 & v88 < 0.9
inc88 <- cong88[, 4]

v90 <- cong90[, 5] / (cong90[, 5] + cong90[, 6])
bad90 <- cong90[, 5] == -9 | cong90[, 6] == -9
v90[bad90] <- NA
contested90 <- v90 > 0.1 & v90 < 0.9
inc90 <- cong90[, 4]

# Plot 7.3
png(paste(PATH, 'plot_ex4_1.png', sep = ''))
v88.hist <- ifelse(v88 < 0.1, 1e-04, ifelse(v88 > 0.9, 1 - 1e-04, v88))
hist(v88.hist, breaks = seq(0, 1, 0.05),
     xlab = "Democratic share of two-party votes", ylab = "", yaxt = "n",
     cex.axis = 1.1, cex.lab = 1.1, cex.main = 1.2, 
     main = "Congress elections (1988)")
dev.off()

# Model fitting
v86.adjusted <- ifelse(v86 < 0.1, 0.25, ifelse(v86 > 0.9, 0.75, v86))
vote.86 <- v86.adjusted[contested88]
incumbency.88 <- inc88[contested88]
vote.88 <- v88[contested88]

install.packages('arm')
library(arm)
fit.88 <- lm(vote.88 ~ vote.86 + incumbency.88)
display(fit.88)

# Figure 7.4
# (a)
par(mfrow = c(1, 1))
png(paste(PATH, 'plot_ex4_2.png', sep = ''))
par(pty = 's', mar = c(5, 5, 4, 1) + 0.1)
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), type = 'n',
     xlab = 'Democratic share for 1986',
     ylab = 'Democratic share for 1988',
     cex.lab = 1)
abline (0, 1, lwd = 0.5)
j.v86 <- ifelse(contested86, v86, jitter(v86, 0.02))
j.v88 <- ifelse(contested88, v88, jitter(v88, 0.02))
points(j.v86[inc88 == 0], j.v88[inc88 == 0], pch = 4)
points(j.v86[inc88 == 1], j.v88[inc88 == 1], pch = 0)
points(j.v86[inc88 == -1], j.v88[inc88 == -1], pch = 15)
mtext("Election results across years", line = 1, cex = 1.2)
dev.off()

# (b)
png(paste(PATH, 'plot_ex4_3.png', sep = ''))
par(pty = 's', mar = c(5, 5, 4, 1) + 0.1)
plot (0, 0, xlim = c(0, 1), ylim = c(0, 1), type = "n",
      xlab = "Democratic share for 1986",
      ylab = "Democratic share for 1988",
      cex.lab = 1)
abline(0, 1, lwd = 0.5)
v86.adjusted <- ifelse(v86 < 0.1, 0.25, ifelse(v86 > 0.9, 0.75, v86))
vote.86 <- v86.adjusted[contested88]
vote.88 <- v88[contested88]
incumbency.88 <- inc88[contested88]
points(vote.86[incumbency.88 == 0], vote.88[incumbency.88 == 0], pch = 4)
points(vote.86[incumbency.88 == 1], vote.88[incumbency.88 == 1], pch = 0)
points(vote.86[incumbency.88 == -1], vote.88[incumbency.88 == -1], pch = 15)
mtext("Adjusted results (correcting uncontested wins to 75%)", line = 1, cex = 1.2)
dev.off()

# Simulation for predicting new data points
incumbency.90 <- inc90
vote.88 <- v88
n.tilde <- length(vote.88)
X.tilde <- cbind(rep(1, n.tilde), vote.88, incumbency.90)

n.sims <- 1000
sim.88 <- sim(fit.88, n.sims)
y.tilde <- array(NA, c(n.sims, n.tilde))
for (s in 1:n.sims) {
  pred <- X.tilde %*% sim.88@coef[s, ]
  ok <- ! is.na(pred)
  y.tilde[s, ok] <- rnorm(sum(ok), pred[ok], sim.88@sigma[s])
}

# Predictive simulation for a nonlinear function of new data
y.tilde.new <- ifelse(y.tilde == 'NaN', 0, y.tilde)
y.tilde.new[is.na(y.tilde)] <- 0
loop <- FALSE
if (loop == FALSE) {
  dems.tilde <- rowSums(y.tilde.new > 0.5)	
} else {
  dems.tilde <- rep (NA, n.sims)
  for (s in 1:n.sims) {
  	dems.tilde[s] <- sum (y.tilde.new[s,] > .5)
  }
}

# Figure 7.5
library(miscTools)
aux <- cbind(sim.88@coef, y.tilde.new, dems.tilde)
aux2 <- aux[1:10, c(1:8, ncol(aux))]
aux3 <- rbind(aux2[1:10, ], colMeans(aux2), colMedians(aux2), apply(aux2, 2, sd))
aux3 <- as.data.frame(aux3)
aux3 <- apply(aux3, 2, round, 4)
aux3[, ncol(aux3)] <- round(aux3[, ncol(aux3)], 1)
rownames(aux3) <- c(paste('simulation', 1:10, sep = ''), 'mean', 'median', 'sd')
colnames(aux3) <- c('beta0', 'beta1', 'beta2', 'pred_y1', 'pred_y2', 'pred_y3',
	                  'pred_y4', 'pred_y5', 'pred_dem_wins')

# Implementation using functions
Pred.88 <- function (X.pred, lm.fit) {
  sim.88 <- sim (lm.fit, 1)
  pred <- X.tilde %*% t(sim.88@coef)
  ok <- ! is.na(pred)
  n.pred <- length(pred)
  y.pred <- rep(NA, n.pred)
  y.pred[ok] <- rnorm (sum(ok), pred[ok], sim.88@sigma)
  return(y.pred)
}

y.tilde <- replicate(1000, Pred.88(X.tilde, fit.88))
dems.tilde <- replicate(1000, Pred.88(X.tilde, fit.88) > 0.5)
# END OF SCRIPT
