#### Sensitivity analysis of a typical RR Lyrae light curve 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

set.seed(3)

source('multiperiod.R')

# Define a typical light curve 
rrlyr <- function(x) { 
  0.34210*cos(x+2.44179)+
    0.17212*cos(2*x+2.48882)+
    0.11765*cos(3*x+2.96421)+
    0.07497*cos(4*x+3.42159)+
    0.04016*cos(5*x+3.96938)+
    0.02059*cos(6*x+4.25067)+
    0.01216*cos(7*x+4.56207)+
    0.00783*cos(8*x+5.27788)
}

# Pick some irregularly spaced points
xs <- seq(0, 4*pi+0.01, 0.01)
noisy_xs <- runif(500, 0, 4*pi)
sigma <- runif(length(noisy_xs), 0.01, 0.1)
noisy_ys <- rrlyr(noisy_xs) + rnorm(length(noisy_xs), 0, sigma)
y_range <- rev(range(rrlyr(xs), noisy_ys) + c(-1,1) * max(sigma))

# Keep originals
noisy_xs_orig <- noisy_xs
noisy_ys_orig <- noisy_ys
sigma_orig <- sigma

# Helper functions for plotting the light curve and its observations
plot_lines <- function() {
  plot(xs/(2*pi), rrlyr(xs), 
       ylim=y_range, axes=FALSE,
       type='l', lty=3, lwd=3, col='darkgray', xaxs='i', xaxt='n', yaxt='n',
       xlab=expression(Phase~phi), ylab=expression(Magnitude~m),
       cex.lab=1.5, tcl=0)
  magaxis(side=1:4, tcl=0.25, labels=c(0,0,0,0), cex.axis=1.25)
}

plot_points <- function() {
  errbar(noisy_xs/(2*pi), noisy_ys, 
         noisy_ys-sigma, noisy_ys+sigma, 
         add=1, col="darkred")
  lines(xs/(2*pi), rrlyr(xs), lty=3, lwd=3, col="darkgray")
}

#setEPS()
#postscript("sensitivity.eps", width=9, height=6, family="Palatino Linotype")
dir.create("plots", showWarnings=FALSE)
cairo_pdf(file.path("plots", "sensitivity.pdf"), 
          width=9, height=8, family=plot_font)
par(mar=c(0, 0, 0, 0), mgp=c(3, 0.25, 0), mfrow=c(3,2), oma=c(5,5,1,1))
for (ii in c(200, 75, 30)) {
  # sample points
  new_is <- sample(1:length(noisy_xs), ii)
  noisy_xs <- noisy_xs_orig[new_is]
  sigma <- sigma_orig[new_is]
  noisy_ys <- noisy_ys_orig[new_is]
  photometry <- data.frame(t=noisy_xs, m=noisy_ys, e=sigma)
  
  # plot least squares fit
  plot_lines(); plot_points()
  text(0.5, y_range[2]+0.15, "Least Squares", pos=1, cex=1.3)
  text(1.5, y_range[2]+0.15, paste(ii, "observations"), pos=1, cex=1.3)
  meven <- fit_single_iterative(photometry, period=2*pi, n_lambda=0)$m_even
  lines(c(phases, 1+phases), rep(meven, 2), lwd=3, col="white")
  lines(c(phases, 1+phases), rep(meven, 2), lwd=2)
  axis(side=2, tcl=0, at=c(-0.5, 0, 0.5))
  if (ii==30) {
    axis(side=1, at=c(0, 0.5, 1, 1.5), labels=c(0, 0.5, 1, 1.5), tcl=0)
    mtext('Magnitude', side=2, line=3, outer=T, cex=1.3)
    mtext('Phase', side=1, line=3, outer=T, cex=1.3)
  }
  
  # plot lasso fit
  plot_lines(); plot_points()
  text(0.5, y_range[2]+0.15, "LASSO", pos=1, cex=1.3)
  text(1.5, y_range[2]+0.15, paste(ii, "observations"), pos=1, cex=1.3)
  meven <- fit_single_iterative(photometry, period=2*pi, n_lambda=1000)$m_even
  lines(c(phases, 1+phases), rep(meven, 2), lwd=4, col="white")
  lines(c(phases, 1+phases), rep(meven, 2), lwd=3)
  if (ii==30) 
    axis(side=1, at=c(0, 0.5, 1, 1.5, 2), labels=c(0, 0.5, 1, 1.5, 2), tcl=0)
}
dev.off()
