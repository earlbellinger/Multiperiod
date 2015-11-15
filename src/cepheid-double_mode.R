#### LASSO fit of a multi-mode Cepheid 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

source('multiperiod.R')

dir.create("plots", showWarnings=FALSE)
cairo_pdf(file.path("plots", "OGLE-SMC-CEP-0408.pdf"), 
          width=9, height=4, family=plot_font)
par(mar=c(5, 5, 1, 1), mgp=c(3, 0.25, 0))

#periods <- c(0.2688471, 0.2173800, 0.1824204)
#fname <- file.path('data', 'OGLE-SMC-CEP-3867.dat')
periods <- c(1.7901765, 1.3140538)
fname <- file.path('..', 'sample_data', 'OGLE-SMC-CEP-0408.dat')
data <- read.table(fname, col.names=c('t', 'm', 'e'))
fit <- fit_multiple(data, periods, n_lambda=1000, Nmax=10)

binwidth <- prod(periods)*10#sum(periods)*10
statsbin <- stats.bin(data$t, data$m, 
                      breaks=seq(data$t[1], data$t[length(data$t)], binwidth))
t_center <- statsbin$centers[which.max(statsbin$stats[1,])]
t_minmax <- c(t_center-binwidth/2, t_center+binwidth/2)
t_range <-  data$t>t_minmax[1] & data$t<t_minmax[2]

new_t <- seq(t_minmax[1]-1, t_minmax[2]+1, 0.1)
predictions <- predict(fit$cvfit, newx=Fourier(new_t, periods, 10)$X)

plot(new_t - data$t[t_range][1], predictions,
     type='l', lty=2, lwd=2, col="#666666",
     xlab="Time [days]",
     ylab="I-band Magnitude",
     xlim=range(data$t[t_range]-data$t[t_range][1]),
     ylim=rev(range(data$m[t_range]) + c(-1, 1) * max(data$e[t_range])),
     tcl=0, cex.lab=1.5)#, xaxt='n', yaxt='n')#, col="darkred")
magaxis(side=1:4, tcl=0.25, labels=c(0,0,0,0))

errbar(data$t[t_range]-data$t[t_range][1], 
       data$m[t_range], 
       data$m[t_range]+data$e[t_range],
       data$m[t_range]-data$e[t_range],
       xaxt='n', yaxt='n', col="darkred", cap=0, add=1)

dev.off()
