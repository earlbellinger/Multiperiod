#### Fit light curves to multiperiodic pulsators 
#### Author: Earl Bellinger ( bellinger@mps.mpg.de ) 
#### Stellar Ages & Galactic Evolution Group 
#### Max-Planck-Institut fur Sonnensystemforschung 

library(glmnet)
library(magicaxis)
library(Hmisc)
library(fields)
library(extrafont)

# Global control parameters
plot_font <- "Palatino"# Linotype"#"FiraSans-Regular"
pdf_plot_width <- 9
pdf_plot_height <- 6
png_plot_width <- 675
png_plot_height <- 450

nlambda <- 100
alpha <- 1
N_max <- 8
phases <- seq(0, 1, 0.001)

# Construct the matrix "X" that we use for linear regression
Fourier <- function(t=phases, p=1, Nmax=N_max) {
  base_freqs <- 1/p
  combos <- expand.grid(rep(list(-Nmax:Nmax), length(base_freqs)))
  freqs <- apply(combos, 1, function(combo) { base_freqs %*% combo })
  freqs <- unique(abs(freqs))
  freqs <- rev(freqs[freqs>0])
  num_freqs <- length(freqs)
  X <- matrix(nrow=length(t), ncol=2*num_freqs) 
  for (ii in 1:num_freqs) {
    X[,2*(ii-1)+1] <- sin(2*pi*t*freqs[ii])
    X[,2*(ii-1)+2] <- cos(2*pi*t*freqs[ii])
  }
  return(list(X=X, freqs=freqs))
}
evenly_spaced <- Fourier()$X # a global variable with the phase defaults 

################################################################################
### Plotting ###################################################################
################################################################################

## Create a pdf or png file in 'save_dir' 
start_devices <- function(save_dir, filename, plot_name, output_fmt) {
  dir.create(save_dir, showWarnings=FALSE)
  id <- sub('.+/', '', sub('\\..+', '', filename))
  if (grepl('pdf', output_fmt)) {
    cairo_pdf(file.path(save_dir, paste0(id, '_', plot_name, '.pdf')), 
              width=pdf_plot_width, height=pdf_plot_height, family=plot_font,
              bg = "transparent")
  } else if (grepl('png', output_fmt)) {
    png(file.path(save_dir, paste0(id, '_', plot_name, '.png')),
        width=png_plot_width, height=png_plot_height, family=plot_font)
  }
  par(cex=1.25)
}

# Shows observed - predicted magnitudes as a function of observed
plot_residuals <- function(m, m_res, err, xlab="Magnitude") {
  errbar(m, m_res, 
         m_res+err/2, 
         m_res-err/2, 
         main='', xaxs='i', tcl=0, cex=0, cap=0,
         xlim=rev(range(m)), ylim=rev(range(m_res)),
         xlab='', ylab=expression(m - hat(m)))
  title(main='Residual plot')
  mtext(xlab, side=1, line=1.5)
  abline(h=0, col='darkred', lty=2)
  magaxis(side=1:4, tcl=0.25, labels=c(0,0,0,0))
}

plot_single <- function(period, photometry, phase, m_even, m_hat, filename) {
  ## Phase plot
  par(las=1, mar=c(2, 4, 3, 1), mgp=c(2.75, 0.25, 0), mfrow=c(3, 1), 
      cex=par()$cex)
  errbar(phase, photometry$m, 
         photometry$m+photometry$e/2, 
         photometry$m-photometry$e/2,
         xlim=c(0,1), ylim=rev(range(photometry$m, m_even)),
         xaxs='i', tcl=0, cex=0, cap=0,
         main='', ylab='Magnitude', xlab="")
  title(main=sub('\\..+', '', filename))
  mtext(paste0("Phase", " (", period, " day period)"), side=1, line=1.5)
  lines(c(phases, 1), c(m_even, m_even[1]), col='darkred')
  magaxis(side=1:4, tcl=0.25, labels=c(0,0,0,0))
  
  ## Residual plot
  m_res <- photometry$m-m_hat
  plot_residuals(photometry$m, m_res, photometry$e)
  
  ## Normal QQ plot
  qqnorm(m_res, pch=3)
  qqline(m_res)
}

plot_multiple <- function(periods, photometry, m_hat, filename, t0s=0) {
  ## Phase plots
  layout(matrix(c(if (length(periods) %% 2 == 0) 1 else c(1,1),
                  2:length(periods), length(periods)+1, length(periods)+2),
                ncol=2, byrow=TRUE))
  par(las=1, mar=c(3, 4, 3, 1), mgp=c(2.75, 0.25, 0), oma=c(0, 0, 0, 0),
      cex=par()$cex)
  for (ii in 1:length(periods)) {
    t0 <- ifelse(length(t0s) == length(periods), t0s[ii], t0)
    phase <- ((photometry$t-t0) %% periods[ii])/periods[ii]
    errbar(phase, photometry$m, 
           photometry$m+photometry$e/2, 
           photometry$m-photometry$e/2,
           xlim=c(0,1), ylim=rev(range(photometry$m)),
           xaxs='i', tcl=0, cex=0, cap=0,
           main='', xlab='', ylab='Magnitude')
    mtext(paste0("Phase", " (", periods[ii], " day period)"), side=1, line=1.5)
    magaxis(side=1:4, tcl=0.25, labels=c(0,0,0,0))
  }
  title(main=sub('\\..+', '', filename), outer=1, line=-2)
  
  ## Residual plot
  m_res <- photometry$m-m_hat
  plot_residuals(photometry$m, m_res, photometry$e)
  
  ## Normal QQ plot
  qqnorm(m_res, pch=3)
  qqline(m_res)
}

################################################################################
### Fitting ####################################################################
################################################################################

fit_single <- function(photometry, period=1, t0=0, n_lambda=nlambda, 
                       alpha.=alpha) {
  phase <- ((photometry$t-t0) %% period)/period
  Fourier_phase <- Fourier(phase, 1)$X
  if (n_lambda > 0) {
    cvfit <- cv.glmnet(Fourier_phase, photometry$m, weights=1/photometry$e, 
                       nlambda=n_lambda, alpha=alpha.)
    m_hat <- predict(cvfit, newx=Fourier_phase, s="lambda.min", exact=TRUE)
    m_even <- predict(cvfit, newx=evenly_spaced, s="lambda.min", exact=TRUE)
    mse <- cvfit$cvm[cvfit$lambda == cvfit$lambda.min]
  } else {
    X <- Fourier_phase
    cvfit <- lm(photometry$m ~ X)
    m_hat <- predict(cvfit)
    X <- evenly_spaced
    m_even <- predict(cvfit, newdata=data.frame(X))
    mse <- sum((m_hat-photometry$m)**2)
  }
  return(list(m_even=m_even, m_hat=m_hat, cvfit=cvfit, phase=phase, mse=mse))
}

fit_single_iterative <- function(photometry, period=1, t0=0, n_lambda=nlambda, 
                                 alpha.=alpha, N_range=1:(N_max*2+1)) {
  phase <- ((photometry$t-t0) %% period)/period
  best <- list(mse=Inf)
  for (Nmax in N_range) {
    if (Nmax > length(photometry$m)) return(best)
    Fourier_phase <- Fourier(phase, 1, Nmax)$X
    
    if (n_lambda > 0) { # Use LASSO
      cvfit <- cv.glmnet(Fourier_phase, photometry$m, weights=1/photometry$e, 
                         nlambda=n_lambda, alpha=alpha.)
      m_hat <- predict(cvfit, newx=Fourier_phase, s="lambda.min", exact=TRUE)
      mse <- cvfit$cvm[cvfit$lambda == cvfit$lambda.min]
      if (mse < best$mse) {
        best <- list(m_even=predict(cvfit, newx=Fourier(phases, 1, Nmax)$X, 
                                    s="lambda.min", exact=TRUE), 
                     m_hat=m_hat, cvfit=cvfit, phase=phase, mse=mse)
      }
      
    } else { # Ordinary least squares with Baart's criterion
      X <- Fourier_phase
      cvfit <- lm(photometry$m ~ X)
      m_hat <- predict(cvfit)
      X <- Fourier(phases, 1, Nmax)$X
      m_even <- predict(cvfit, newdata=data.frame(X))
      mse <- sum((m_hat-photometry$m)**2)
      best <- list(m_even=m_even, m_hat=m_hat, cvfit=cvfit, 
                   phase=phase, mse=mse)
      resids <- as.numeric((m_hat - photometry$m)[order(phase)])
      mu <- mean(resids)
      nums <- 0
      dens <- 0
      for (ii in 1:(length(photometry$m)-1)) {
        nums <- nums + (resids[ii]-mu)*(resids[ii+1]-mu)
        dens <- dens + (resids[ii]-mu)**2
      }
      cutoff <- (2*(length(photometry$m)-1))**(-1/2)
      if (abs(nums/dens) <= cutoff) return(best)
    }
  }
  return(best)
}

fit_multiple <- function(photometry, periods, n_lambda=nlambda, alpha.=alpha,
                         Nmax=N_max) {
  Fourier_space <- Fourier(photometry$t, periods, Nmax)
  cvfit <- cv.glmnet(Fourier_space$X, photometry$m, weights=1/photometry$e, 
                     nlambda=n_lambda, alpha=alpha.)
  m_hat <- predict(cvfit, newx=Fourier_space$X, s="lambda.min", exact=TRUE)
  mse <- cvfit$cvm[cvfit$lambda == cvfit$lambda.min]
  return(list(m_hat=m_hat, cvfit=cvfit, freqs=Fourier_space$freqs,
              coefs=coef(cvfit), mse=mse))
}

obj_f <- function(photometry, periods, 
                  lambda=NULL, n_lambda=nlambda, alpha.=alpha) {
  Fourier_space <- Fourier(photometry$t, periods)
  cvfit <- cv.glmnet(Fourier_space$X, photometry$m, weights=1/photometry$e, 
                     nlambda=n_lambda, alpha=alpha.)
  coefs <- coef(cvfit, s = "lambda.min")
  sum((photometry$m - coefs[1] - (Fourier_space$X %*% coefs[-1]))^2) + 
    cvfit$lambda.min * ifelse(alpha.==1, sum(abs(coefs[-1])), sum(coefs[-1]^2))
}

################################################################################
### Interface ##################################################################
################################################################################

fit_lightcurve <- function(filename, periods=1, t0s=0, show_plot=TRUE, 
                           save_dir=NULL, output_fmt='.png', n_lambda=nlambda) {
  if (!file.exists(filename)) 
    return("Cannot find file")
  photometry <- read.table(filename, col.names=c('t', 'm', 'e'))
  
  if (show_plot && !is.null(save_dir))
    start_devices(save_dir, filename, 'phase', output_fmt)
  
  if (length(periods) == 1) {
    fit <- fit_single(photometry, periods, t0=t0s[1], n_lambda)
    if (show_plot) 
      plot_single(periods, photometry, fit$phase, 
                  fit$m_even, fit$m_hat, filename)
  } else { 
    fit <- fit_multiple(photometry, periods, n_lambda)
    if (show_plot) 
      plot_multiple(periods, photometry, fit$m_hat, filename, t0s)
  }
  
  if (show_plot && !is.null(save_dir)) dev.off()
  return(fit)
}

test <- function() {
  fit_lightcurve('data/OGLE-LMC-CEP-0002.dat', 3.118120, 2171.239,
                 save_dir='multiplots', output_fmt='.png')
  
  fit_lightcurve('data/OGLE-SMC-CEP-0408.dat', 
                 c(1.7901765, 1.3140538), c(624.75372, 624.97043),
                 save_dir='multiplots', output_fmt='.png')
  
  fit_lightcurve('data/OGLE-SMC-CEP-3867.dat', 
                 c(0.2688471, 0.2173800, 0.1824204), 
                 c(2104.60491, 2104.81172, 2104.72340),
                 save_dir='multiplots', output_fmt='.png')
}
