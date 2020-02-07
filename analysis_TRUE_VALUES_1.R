#------------------------------------------------------
#------------------------------------------------------
# analysis_TRUE_VALUES_1.R
# - Obtains (estimated) true values of the cumulative coefficients of the correctly-specified MSM, and of the causal estimands of interest by fitting the MSM to the large trial data generated in dat_sim_TRUE_VALUES_1.R. 
#------------------------------------------------------
#------------------------------------------------------

#--------------------------
#packages
#--------------------------

library(survival)
library(timereg)

#--------------------------
#times at which we obtain estimated cumulative coefficients of the MSM and estimated survival probabilities
#--------------------------

t.hor=seq(0,5,0.01)

#--------------------------
#analysis of n.sim large trials
#--------------------------

n.sim.true=1000

#----
#store results

coef_0=matrix(nrow=n.sim.true,ncol=length(t.hor))
coef_A=matrix(nrow=n.sim.true,ncol=length(t.hor))
coef_Alag1=matrix(nrow=n.sim.true,ncol=length(t.hor))
coef_Alag2=matrix(nrow=n.sim.true,ncol=length(t.hor))
coef_Alag3=matrix(nrow=n.sim.true,ncol=length(t.hor))
coef_Alag4=matrix(nrow=n.sim.true,ncol=length(t.hor))

surv1=matrix(nrow=n.sim.true,ncol=length(t.hor))
surv0=matrix(nrow=n.sim.true,ncol=length(t.hor))

#----
#start of simulation loop

for(i in 1:n.sim.true){
  
  print(i)

  #-----
  #simulate data using dat_sim_TRUE_VALUES_1
  
  source("dat_sim_TRUE_VALUES_1.R")
  
  #-----
  #Fit the correctly specified MSM 
  #No inverse probability weights are needed because there is no time-dependent confounding here
  #Note that the coefficient for Alag1 is 0 for t<=1; the coefficient for Alag2 is 0 for t<=2, etc.
  #The aalen function doesn't like this, and if we try and fit a single model for 0<t<=5 it only gives estimates for 4<t<=5, where all coefficients are non-zero
  #Hence we fit the model separately in a series of time periods: 0<t<=1, 1<t<=2,..., 4<t<=5
  
  ah.p0=aalen(Surv(time,time.stop,event)~A,data=dat.long[dat.long$time==0,],n.sim=0)
  ah.p1=aalen(Surv(time,time.stop,event)~A+Alag1,data=dat.long[dat.long$time==1,],n.sim=0)
  ah.p2=aalen(Surv(time,time.stop,event)~A+Alag1+Alag2,data=dat.long[dat.long$time==2,],n.sim=0)
  ah.p3=aalen(Surv(time,time.stop,event)~A+Alag1+Alag2+Alag3,data=dat.long[dat.long$time==3,],n.sim=0)
  ah.p4=aalen(Surv(time,time.stop,event)~A+Alag1+Alag2+Alag3+Alag4,data=dat.long[dat.long$time==4,],n.sim=0)
  
  #-----
  #Create a step function for each estimated cumulative coefficient, enabling us to obtain the estimated cumulative coefficient at any time 
  
  maxrow.p0=dim(ah.p0$cum)[1]
  maxrow.p1=dim(ah.p1$cum)[1]
  maxrow.p2=dim(ah.p2$cum)[1]
  maxrow.p3=dim(ah.p3$cum)[1]
  maxrow.p4=dim(ah.p4$cum)[1]
  
  ah.stepfun.0=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                        c(0,ah.p0$cum[,2],
                          ah.p0$cum[maxrow.p0,2]+ah.p1$cum[-1,2],
                          ah.p0$cum[maxrow.p0,2]+ah.p1$cum[maxrow.p1,2]+ah.p2$cum[-1,2],
                          ah.p0$cum[maxrow.p0,2]+ah.p1$cum[maxrow.p1,2]+ah.p2$cum[maxrow.p2,2]+ah.p3$cum[-1,2],
                          ah.p0$cum[maxrow.p0,2]+ah.p1$cum[maxrow.p1,2]+ah.p2$cum[maxrow.p2,2]+ah.p3$cum[maxrow.p3,2]+ah.p4$cum[-1,2]))
  
  ah.stepfun.A=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                        c(0,ah.p0$cum[,3],
                          ah.p0$cum[maxrow.p0,3]+ah.p1$cum[-1,3],
                          ah.p0$cum[maxrow.p0,3]+ah.p1$cum[maxrow.p1,3]+ah.p2$cum[-1,3],
                          ah.p0$cum[maxrow.p0,3]+ah.p1$cum[maxrow.p1,3]+ah.p2$cum[maxrow.p2,3]+ah.p3$cum[-1,3],
                          ah.p0$cum[maxrow.p0,3]+ah.p1$cum[maxrow.p1,3]+ah.p2$cum[maxrow.p2,3]+ah.p3$cum[maxrow.p3,3]+ah.p4$cum[-1,3]))
  
  ah.stepfun.Alag1=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                            c(0,rep(0,dim(ah.p0$cum)[1]),
                              ah.p1$cum[-1,4],
                              ah.p1$cum[maxrow.p1,4]+ah.p2$cum[-1,4],
                              ah.p1$cum[maxrow.p1,4]+ah.p2$cum[maxrow.p2,4]+ah.p3$cum[-1,4],
                              ah.p1$cum[maxrow.p1,4]+ah.p2$cum[maxrow.p2,4]+ah.p3$cum[maxrow.p3,4]+ah.p4$cum[-1,4]))
  
  ah.stepfun.Alag2=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                            c(0,rep(0,dim(ah.p0$cum)[1]),
                              rep(0,dim(ah.p1$cum)[1]-1),
                              ah.p2$cum[-1,5],
                              ah.p2$cum[maxrow.p2,5]+ah.p3$cum[-1,5],
                              ah.p2$cum[maxrow.p2,5]+ah.p3$cum[maxrow.p3,5]+ah.p4$cum[-1,5]))
  
  ah.stepfun.Alag3=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                            c(0,rep(0,dim(ah.p0$cum)[1]),
                              rep(0,dim(ah.p1$cum)[1]-1),
                              rep(0,dim(ah.p2$cum)[1]-1),
                              ah.p3$cum[-1,6],
                              ah.p3$cum[maxrow.p3,6]+ah.p4$cum[-1,6]))
  
  ah.stepfun.Alag4=stepfun(c(ah.p0$cum[,1], ah.p1$cum[-1,1],ah.p2$cum[-1,1],ah.p3$cum[-1,1],ah.p4$cum[-1,1]),
                            c(0,rep(0,dim(ah.p0$cum)[1]),
                              rep(0,dim(ah.p1$cum)[1]-1),
                              rep(0,dim(ah.p2$cum)[1]-1),
                              rep(0,dim(ah.p3$cum)[1]-1),
                              ah.p4$cum[-1,7]))
  
  #----
  #estimated cumulative coefficients from the MSM at times t.hor
  
  coef_0[i,]=ah.stepfun.0(t.hor)
  coef_A[i,]=ah.stepfun.A(t.hor)
  coef_Alag1[i,]=ah.stepfun.Alag1(t.hor)
  coef_Alag2[i,]=ah.stepfun.Alag2(t.hor)
  coef_Alag3[i,]=ah.stepfun.Alag3(t.hor)
  coef_Alag4[i,]=ah.stepfun.Alag4(t.hor)
  
  #----
  #estimated survival probabilities at times t.hor in the 'always treated' (surv1) and in the 'never treated' (surv0)
  
  surv1[i,]=exp(-(ah.stepfun.0(t.hor)+ah.stepfun.A(t.hor)+ah.stepfun.Alag1(t.hor)+ah.stepfun.Alag2(t.hor)+ah.stepfun.Alag3(t.hor)+ah.stepfun.Alag4(t.hor)))
  surv0[i,]=exp(-(ah.stepfun.0(t.hor)))
}

#end of simulation loop
#----

#--------------------------
#True values of cumulative coefficients and survival probabilities
#Taken to be the average of the n.sim.true estimates from the large trials
#--------------------------

coef_0_true=colMeans(coef_0)
coef_A_true=colMeans(coef_A)
coef_Alag1_true=colMeans(coef_Alag1)
coef_Alag2_true=colMeans(coef_Alag2)
coef_Alag3_true=colMeans(coef_Alag3)
coef_Alag4_true=colMeans(coef_Alag4)

surv1_true=colMeans(surv1)
surv0_true=colMeans(surv0)

