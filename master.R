#------------------------------------------------------
#------------------------------------------------------
# master.R
# - Implements the simulation study as described in Section 6 of the paper.
# - Summarises the simulation results to create the numbers shown in Tables 1 and 2 and the plots in Figures 2 and 3 of the paper.
# - Saves in the results folder: table1.csv,table2.csv,
# figure2_1.png, figure2_2.png, figure2_3.png, figure2_4.png, figure2_5.png, figure2_6.png, figure 3.
# - Figure 2 is saved as 6 subplots.
# Created by Ruth Keogh
#------------------------------------------------------
#------------------------------------------------------

#--------------------------
#packages
#--------------------------

library(survival)
library(timereg)
library(ggplot2)
library(gridExtra)
library(tidyverse)

#--------------------------
# functions.R
#--------------------------
# - Loads the functions needed.

source("functions.R")

#--------------------------
# analysis.R
#--------------------------
# - Simulates longitudinal and time-to-event data using the algorithm described in Section 5 of the paper.
# - Data are simulated using the file dat_sim.R
# - Fits a correctly specified MSM (using IPTW) to the simulated data. The MSM is an additive hazard model. See Section 6.1 of the paper.
# - Obtains estimates of cumulative coefficients from the MSM. 
# - Obtains estimates of survival probabilities under the treatment regimes 'always treated' and 'never treated'. 
# - Repeats over 1000 simulated data sets. 
# - can be replaced by analysis_CENS.R

source("analysis.R")

#--------------------------
# analysis_TRUE_VALUES_1.R
#--------------------------
# - Obtains (estimated) true values of the cumulative coefficients of the correctly-specified MSM, 
# and of the causal estimands of interest by fitting the MSM to the large trial data generated in dat_sim_TRUE_VALUES_1.R. 
# - Can be replaced by analysis_TRUE_VALUES_2.R

source("analysis_TRUE_VALUES_1.R")

#--------------------------
# sim_results.R
#--------------------------
# - Summarises the simulation results to create Tables 1 and 2 and Figures 2 and 3 in the paper.
# - Saves in the results folder: table1.csv,table2.csv,
# figure2_1.png, figure2_2.png, figure2_3.png, figure2_4.png, figure2_5.png, figure2_6.png, figure 3.

source("sim_results.R")
