

draw <- function(m, pic.file) {
  alpha <- seq(0, 1, 1/10000)
  f1 <- 1 - (1 - alpha) ** (1 / m)
  f2 <- alpha / m
  f3 <- alpha
  png(paste('~/Desktop/bgse/projects/stats/', pic.file, sep = ''))
  plot(0, xlim = c(0,1), ylim = c(0,1),
       main = paste('Function shape (m = ', m, ')', sep = ''),
       xlab = 'Input', ylab = 'Output', pch = '')
  lines(alpha, f1, col = 'blue', lwd = 2)
  lines(alpha, f2, col = 'red', lwd = 2)
  lines(alpha, f3, col = 'darkgreen', lwd = 2)
  legend('topleft', c('f_1', 'f_2', 'f_3'),
                    col = c('blue','red', 'darkgreen'), lty = 1, lwd = 2)
  dev.off()
}

draw(m = 1, 'pic1.png')
draw(m = 2, 'pic2.png')
draw(m = 5, 'pic3.png')
draw(m = 10, 'pic4.png')
draw(m = 20, 'pic5.png')



