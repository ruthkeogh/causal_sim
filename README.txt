#######################################
#######################################
This project contains R code for implementation of the simulation study described in the paper: 

TITLE: Simulating longitudinal data from marginal structural models using the additive hazards model. 
JOURNAL: Biometrical Journal 2021. 
AUTHORS: Ruth Keogh, Shaun Seaman SR, Jon Michael Gran, Stijn Vansteelandt. 

Ruth Keogh was responsible for writing the code. 
Please address questions, comments, and reporting of bugs, to Ruth Keogh: ruth.keogh@lshtm.ac.uk

R version 3.6.3 (2020-02-29)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19042)
#######################################
#######################################


The following code files are provided and generate Tables 1 and 2 and Figures 2 and 3 from the paper, as also shown in simulation_replication.pdf. It is recommended to start by reading the short document: simulation_replication.pdf.

1. 
master.R
- Runs functions.R, analysis.R, analysis_TRUE_VALUES_1.R, sim_results.R.
- Outputs Tables 1 and 2 to the 'results' folder (table1.csv,table2.csv).
- Outputs Figures 2 and 3 to the 'results' folder (figure2_1.png, figure2_2.png, figure2_3.png, figure2_4.png figure2_5.png, figure2_6.png, figure 3). Figure 2 is saved as 6 subplots. 

2. 
functions.R
- defines several functions used in the analysis and to produce the results summary. 

3.
dat_sim.R
- Simulates longitudinal and time-to-event data using the algorithm described in Section 5 of the paper.
- Is called in analysis.R.

4.
analysis.R
- Calls dat_sim.R
- Fits a correctly specified MSM (using IPTW) to the simulated data. The MSM is an additive hazard model. See Section 6.1 of the paper.
- Obtains estimates of cumulative coefficients from the MSM. 
- Obtains estimates of survival probabilities under the treatment regimes 'always treated' and 'never treated'. 
- Repeats over 1000 simulated data sets. 

5.
dat_sim_TRUE_VALUES_1.R
- Generates data from a large trial to obtain true values of the cumulative coefficients of the correctly-specified MSM, and of the causal estimands of interest (survival probabilities under the treatment regimes 'always treated' and 'never treated'). 
- See Section 6.2 of the paper.
- Is called in analysis_TRUE_VALUES_1.R.

6.
analysis_TRUE_VALUES_1.R
- Calls dat_sim_TRUE_VALUES_1.R
- Obtains (estimated) true values of the cumulative coefficients of the correctly-specified MSM, and of the causal estimands of interest by fitting the MSM to the large trial data generated in dat_sim_TRUE_VALUES_1.R.
- See Section 6.2 of the paper

7.
sim_results.R
- Summarises the simulation results to create Tables 1 and 2 and Figures 2 and 3 of the paper.
- The tables and figures are saved to the results folder. 

#######################################
Files 8-11 are extra and include extensions of the above files to other scenarios discussed in the paper.
#######################################

8. dat_sim_TRUE_VALUES_2.R
- An alternative to dat_sim_TRUE_VALUES_1.R
- Unlike dat_sim_TRUE_VALUES_1.R, this file only simulates data under two treatment regimes: 'always treated' and 'never treated'.
- This saves computing time. 
- It is accompanied by analysis_TRUE_VALUES_2.R

9. analysis_TRUE_VALUES_2.R
- Can replace analysis_TRUE_VALUES_1.R in master.R.
- Obtains (estimated) true values of survival probabilities under the treatment regimes 'always treated' and 'never treated' using the large trial data generated in dat_sim_TRUE_VALUES_2.R.  
- Survival probabilities are directly estimated using simple proportions, since there is only administrative censoring in the data.
- This is just an another way of obtaining the true values for the survival probabilities. It is faster than using dat_sim_TRUE_VALUES_1.R and analysis_TRUE_VALUES_1.R because it doesn't involve fitting a model. See Section 6.2. of the paper. 
- Unlike analysis_TRUE_VALUES_1.R, this file does not obtain (estimated) true values of the cumulative coefficients. Therefore it can be used to create Table 2 and Figure 3, but not Table 1 or Figure 2.

10. dat_sim_CENS.R
- Can replace dat_sim.R in analysis.R.
- Simulates longitudinal and time-to-event data using the algorithm described in Section 5, with additional right-censoring that depends on the most recent values of L (it could also depend on A).
- The parameter settings reult in around 30% of individuals being censored.

11. dat_sim_GEN.R
- Can replace dat_sim.R in analysis.R.
- An example of simulating longitudinal and time-to-event data using the algorithm described in the Supplementary Materials Section A5, accommodating a more complex (non- exponential) form for the conditional additive hazard model. 



