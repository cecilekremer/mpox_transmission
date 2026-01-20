# System requirements
The code used in this manuscript was developed using R version 4.3.2.
To install R, we refer to https://cran.r-project.org/doc/manuals/r-release/R-admin.html

# R code
The following R scripts are included in this repository:
* 00_descriptives: descriptive analyses of combined data (Goma + Kamituga)
* 01a_serial_interval_observed.R: to obtain parameter estimates for the observed serial intervals
* 01b_serial_interval_regMCMC.R: to estimate parameters of the Bayesian linear regression model for the serial interval
* fun_reg_si.R: function to perform MCMC for the regression model
* fun_network.R: function to sample transmission trees
* 02_incubation_period.R: to estimate parameters of the incubation period + Bayesian linear regression model for the incubation period
* stan_incub_reg.stan: Stan model for the Bayesian linear regression model for the incubation period
* 03_ct_analysis.R: script to perform the analysis of Ct data
