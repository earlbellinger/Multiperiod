# multiperiod
R and MATLAB packages for the multiperiodic Fourier decomposition. 

The file `src/multiperiod.R` is an R implementation of the multiperiodic least squares and LASSO Fourier decompositions, which are useful for reconstructing periodic signals from noisy observations that have been sampled irregularly in time. 

The file `src/sensitivity.R` performs a sensitivity comparison between the least squares and LASSO fitting methods by applying them each to a simulated RR Lyrae light curve with fewer and fewer data points. 

The file `src/simulation.m` is a (very minimal) MATLAB implementation of the least squares and LASSO Fourier fitting methods. It generates a simulated multiperiodic oscillator and fits it with both methods. 

Why implementations in multiple languages? Well, to minimize the chance that we've done something wrong, of course! 

--

This repository hosts the source code for the figures and TeX of 

Bellinger, E. P., Wysocki, D., Kanbur, S. M. (2015). Measuring amplitudes of harmonics and combination frequencies in variable stars. *Communications from the Konkoly Observatory of the Hungarian Academy of Sciences*, 105.

If you use this source code in research that leads to publication, please consider citing that paper. 
