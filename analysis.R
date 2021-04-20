#------------------------------------------------------
#------------------------------------------------------
# analysis.R
# - Fits a correctly specified MSM (using IPTW) to the simulated data. The MSM is an additive hazard model. See Section 6.1 of the paper.
# - Obtains estimates of cumulative coefficients from the MSM. 
# - Obtains estimates of survival probabilities under the treatment regimes 'always treated' and 'never treated'. 
# - Repeats over 1000 simulated data sets. 
# Created by Ruth Keogh
#------------------------------------------------------
#------------------------------------------------------

#----
#set seed
set.seed(443792)

#--------------------------
#times at which we obtain estimated cumulative coefficients of the MSM and estimated survival probabilities
#--------------------------

t.hor=seq(0,5,0.01)

#--------------------------
#analysis in n.sim simulated data sets
#--------------------------

n.sim=1000

#----
#store results

coef_0=matrix(nrow=n.sim,ncol=length(t.hor))
coef_A=matrix(nrow=n.sim,ncol=length(t.hor))
coef_Alag1=matrix(nrow=n.sim,ncol=length(t.hor))
coef_Alag2=matrix(nrow=n.sim,ncol=length(t.hor))
coef_Alag3=matrix(nrow=n.sim,ncol=length(t.hor))
coef_Alag4=matrix(nrow=n.sim,ncol=length(t.hor))

surv1=matrix(nrow=n.sim,ncol=length(t.hor))
surv0=matrix(nrow=n.sim,ncol=length(t.hor))

#----
#start of simulation loop

for(i in 1:n.sim){
  
  print(i)
  
  #-----
  #simulate data using dat_sim
  #could be replaced with dat_sim_GEN.R
  
  source("dat_sim.R")
  
  #------------------
  #stabilized weights
  #------------------
  
  wt.mod=glm(A~L+Alag1,family="binomial",data=dat.long)
  pred.wt=predict(wt.mod,type = "response")
  
  dat.long$wt=ifelse(dat.long$A==1,pred.wt,1-pred.wt)
  dat.long$wt.cum=ave(dat.long$wt,dat.long$id,FUN=cumprod)
  
  wt.mod.num=glm(A~Alag1,family="binomial",data=dat.long)
  pred.wt.num=predict(wt.mod.num,type = "response")
  
  dat.long$wt.num=ifelse(dat.long$A==1,pred.wt.num,1-pred.wt.num)
  dat.long$wt.cum.num=ave(dat.long$wt.num,dat.long$id,FUN=cumprod)
  
  dat.long$ipw.s=dat.long$wt.cum.num/dat.long$wt.cum
  
  #-----
  #Fit the correctly specified MSM using IPTW
  #See notes in analysis_TRUE_VALUES.R which explain why we fit the model separately by time period
  
  ah.p0=aalen(Surv(time,time.stop,event)~A,data=dat.long[dat.long$time==0,],n.sim=0,weights = dat.long[dat.long$time==0,]$ipw.s)
  ah.p1=aalen(Surv(time,time.stop,event)~A+Alag1,data=dat.long[dat.long$time==1,],n.sim=0,weights = dat.long[dat.long$time==1,]$ipw.s)
  ah.p2=aalen(Surv(time,time.stop,event)~A+Alag1+Alag2,data=dat.long[dat.long$time==2,],n.sim=0,weights = dat.long[dat.long$time==2,]$ipw.s)
  ah.p3=aalen(Surv(time,time.stop,event)~A+Alag1+Alag2+Alag3,data=dat.long[dat.long$time==3,],n.sim=0,weights = dat.long[dat.long$time==3,]$ipw.s)
  ah.p4=aalen(Surv(time,time.stop,event)~A+Alag1+Alag2+Alag3+Alag4,data=dat.long[dat.long$time==4,],n.sim=0,weights = dat.long[dat.long$time==4,]$ipw.s)
  
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
