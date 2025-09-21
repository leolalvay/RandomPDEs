This folder solves the item 4 of section 1.4 Numerical Tasks

The main function runs Importance sampling Monte Carlo and uses mu and 
Sigma which are obtained like this:

1. optimize_mu.m 
Solves the optimization problem for \hat{\mu} in the overleaf documen. The
result mu is a vector indicating where the function f attains its maximun, 

2. plot_f_vs_gaussian.m
Plots the contours of the function f agains the contours of the bivariate normal
with mean mu and covariance matrix Sigma (Figure 5 in the overleaf). 
The mu values is taken from the  optimize_mu.m function. The Sigma values 
are precisely tunned using this plot_f_vs_gaussian.m function.

The function: functionfplot.m plots the Figure 4 in the overleaf document.
This is just for analysis.

The outpust of main function are the convergence rates comparison between
Monte Carlo and Importance Sampling (Figure 6 in the Overleaf). These figures
are inside folder Exercise14 with names: convergence_K3 and convergence_K6