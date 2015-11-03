
PATH <- '/Users/miquel/Desktop/bgse/courses/term1/smi/'
DATADIR <- paste(PATH, 'datasets/', sep = '')
load(file = paste(DATADIR, 'labelled_suspicious.RData', sep = ''))

t <- labelled_suspicious[, 't']
phi1 <- labelled_suspicious[, 'phi1']
phi2 <- labelled_suspicious[, 'phi2']

x <- abs((phi1 - mean(phi1)) / sd(phi1)) < 2
plot(t[x], phi1[x])
lines(seq(0, 1, 1), c(mean(phi1[t == 0]), mean(phi1[t == 1])))
lines(seq(0, 1, 1), c(median(phi1[t == 0]), median(phi1[t == 1])))

plot(t, phi2)

###
cond <- function(A) {
  ev <- eigen(A)[['values']]
  return(max(ev) * min(ev) ** (-1))
}

