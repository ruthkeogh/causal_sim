Files for implementation of the simulated methods described in the paper: 
Keogh RH, Seaman SR, Gran JM, Vansteelandt S. Simulating longitudinal data from marginal structural models using the additive hazards model. 2020.

dat_sim.R
- Simulates longitudinal and time-to-event data using the algorithm described in Section 5. 

analysis.R
- Fits a correctly specified MSM (using IPTW) to the simulated data. The MSM is an additive hazard model. See Section 6.1 of the paper.
- Obtains estimates of cumulative coefficients from the MSM. 
- Obtains estimates of survival probabilities under the treatment regimes 'always treated' and 'never treated'. 
- Repeats over 1000 simulated data sets. 

dat_sim_TRUE_VALUES_1.R
- Generates data from a large trial to obtain true values of the cumulative coefficients of the correctly-specified MSM, and of the causal estimands of interest (survival probabilities under the treatment regimes 'always treated' and 'never treated'). See Section 6.2. of the paper

analysis_TRUE_VALUES_1.R
- Obtains (estimated) true values of the cumulative coefficients of the correctly-specified MSM, and of the causal estimands of interest by fitting the MSM to the large trial data generated in dat_sim_TRUE_VALUES_1.R. 

sim_results.R
- Summarises the simulation results to create the numbers shown in Tables 1 and 2 and the plots in Figures 2 and 3 of the paper.

dat_sim_TRUE_VALUES_2.R
- Generates data from a large trial focusing on only two treatment regimes: 'always treated' and 'never treated'.

analysis_TRUE_VALUES_2.R
- Obtains (estimated) true values of survival probabilities under the treatment regimes 'always treated' and 'never treated' using the large trial data generated in dat_sim_TRUE_VALUES_2.R.  
- Survival probabilities are directly estimated using simple proportions, since there is only administrative censoring in the data.
- This is just an another way of obtaining the true values for the survival probabilities. It is faster than using dat_sim_TRUE_VALUES_1.R and analysis_TRUE_VALUES_1.R because it doesn't involve fitting a model. See Section 6.2. of the paper. 

dat_sim_GEN.R
- An example of simulating longitudinal and time-to-event data using the algorithm described in the Supplementary Materials Section A5, accommodating a more complex (non-exponential) form for the conditional additive hazard model. 
