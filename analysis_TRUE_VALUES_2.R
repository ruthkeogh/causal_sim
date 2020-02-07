#------------------------------------------------------
#------------------------------------------------------
# analysis_TRUE_VALUES_2.R
# - Obtains (estimated) true values of survival probabilities under the treatment regimes 'always treated' and 'never treated' using the large trial data generated in dat_sim_TRUE_VALUES_2.R.  
# - Survival probabilities are directly estimated using simple proportions, since there is only administrative censoring in the data.
# - This is just an another way of obtaining the true values for the survival probabilities. It is faster than using dat_sim_TRUE_VALUES_1.R and analysis_TRUE_VALUES_1.R because it doesn't involve fitting a model. See Section 6.2. of the paper. 
#------------------------------------------------------
#------------------------------------------------------


#--------------------------
#times at which we obtain estimated cumulative coefficients of the MSM and estimated survival probabilities
#--------------------------

t.hor=seq(0,5,0.01)

t.hor=seq(0,5,0.01)#time horizon for survival probabilities

#--------------------------
#analysis of n.sim large trials
#--------------------------

n.sim.true=1000

#----
#store results

surv1=matrix(nrow=n.sim.true,ncol=length(t.hor))
surv0=matrix(nrow=n.sim.true,ncol=length(t.hor))

#----
#start of simulation loop

for(i in 1:n.sim.true){
  
  print(i)
  
  #-----
  #simulate data using dat_sim_TRUE_VALUES_2
  
  source("dat_sim_TRUE_VALUES_2.R")

  #----
  #estimated survival probabilities at times t.hor in the 'always treated' (surv1) and in the 'never treated' (surv0)

  surv1[i,]=sapply(1:length(t.hor),function(x){mean(dat.1$T.obs>t.hor[x])})
  
  surv0[i,]=sapply(1:length(t.hor),function(x){mean(dat.0$T.obs>t.hor[x])})
  
}

#end of simulation loop
#----

#--------------------------
#True values of survival probabilities
#Taken to be the average of the n.sim.true estimates from the large trials
#--------------------------

surv1_true=colMeans(surv1)
surv0_true=colMeans(surv0)
