#------------------------------------------------------
#------------------------------------------------------
# dat_sim_CENS.R
# - Simulates longitudinal and time-to-event data using the algorithm described in Section 5. 
# - with additional right-censoring that depends on the most recent values of L (it could also depend on A)
# - the parameter settings here reult in around 30% of individuals being censored
# Created by Ruth Keogh
#------------------------------------------------------
#------------------------------------------------------

#----
#sample size 

n=5000

#----
#number of visits (K+1)

n.visit=5

#-------------------
#data generation
#-------------------

#----
#first generate the time-dependent A and L

A=matrix(nrow=n,ncol=n.visit)
L=matrix(nrow=n,ncol=n.visit)

U=rnorm(n,0,0.1)

L[,1]=rnorm(n,0+U,1)
A[,1]=rbinom(n,1,expit(-2+0.5*L[,1]))
for(k in 2:n.visit){
  L[,k]=rnorm(n,0.8*L[,k-1]-A[,k-1]+0.1*(k-1)+U,1)
  A[,k]=rbinom(n,1,expit(-2+0.5*L[,k]+A[,k-1]))
}

#----
#generate event times T.obs, and event indicators D.obs

T.obs=rep(NA,n)

for(v in 1:n.visit){
  
  u.t=runif(n,0,1)
  haz=0.7-0.2*A[,v]+0.05*L[,v]+0.05*U
  new.t=-log(u.t)/haz
  T.obs=ifelse(is.na(T.obs) & new.t<1 & haz>0,v-1+new.t,T.obs)#the haz>0 is just used to deal with tiny possibility (under this data generating mechanism) the hazard could go negative. 
  
}

#----
#generate additional right-censoring
#people can be censored at any visit time
#an alternative would be to generate censoring times in continuous time
#C=0 if a person is observed at a given visit and 1 if they are not observed (i.e. are censored at that visit)

C=matrix(nrow=n,ncol=n.visit)
C[,1]=1
for(k in 2:n.visit){
  C[,k]=rbinom(n,1,expit(0+1*L[,k]))
}

T.cens=5 #censoring time is the earliest time at which censoring occurs
T.cens=ifelse(C[,2]==1,1,T.cens)#censored at time 1
T.cens=ifelse(C[,2]==0 & C[,3]==1,2,T.cens)#censored at time 2
T.cens=ifelse(C[,2]==0 & C[,3]==0 & C[,4]==1,3,T.cens)#censored at time 3
T.cens=ifelse(C[,2]==0 & C[,3]==0 & C[,4]==0 & C[,5]==1,4,T.cens)#censored at time 4

D.obs=ifelse(is.na(T.obs),0,1)#new censoring indicator
D.obs=ifelse(T.cens<T.obs & !is.na(T.obs),0,D.obs)
T.obs=ifelse(D.obs==0,T.cens,T.obs)

# table(D.obs)
# table(T.obs[D.obs==0])

#----
#The above data are in 'wide' format (1 row per individual). Reshape into 'long' format (multiple rows per individual: 1 row for each visit)

L.dat=as.data.frame(L)
names(L.dat)=paste0("L.",0:4)

A.dat=as.data.frame(A)
names(A.dat)=paste0("A.",0:4)

Alag1.dat=as.data.frame(cbind(rep(0,n),A.dat[,1:4]))
Alag2.dat=as.data.frame(cbind(rep(0,n),rep(0,n),A.dat[,1:3]))
Alag3.dat=as.data.frame(cbind(rep(0,n),rep(0,n),rep(0,n),A.dat[,1:2]))
Alag4.dat=as.data.frame(cbind(rep(0,n),rep(0,n),rep(0,n),rep(0,n),A.dat[,1]))

names(Alag1.dat)=paste0("Alag1.",0:4)
names(Alag2.dat)=paste0("Alag2.",0:4)
names(Alag3.dat)=paste0("Alag3.",0:4)
names(Alag4.dat)=paste0("Alag4.",0:4)

dat=data.frame(id=1:n,T.obs,D.obs,A.dat,Alag1.dat,Alag2.dat,Alag3.dat,Alag4.dat,L.dat,U)

dat.long=reshape(data = dat,varying=c(paste0("A.",0:4),paste0("Alag1.",0:4),paste0("Alag2.",0:4),paste0("Alag3.",0:4),paste0("Alag4.",0:4),paste0("L.",0:4)),direction="long",idvar="id")

dat.long=dat.long[order(dat.long$id,dat.long$time),]

dat.long$time.stop=dat.long$time+1

dat.long=dat.long[dat.long$time<dat.long$T.obs,]

dat.long$time.stop=ifelse(dat.long$time.stop>dat.long$T.obs,dat.long$T.obs,dat.long$time.stop)

dat.long$event=ifelse(dat.long$time.stop==dat.long$T.obs & dat.long$D.obs==1,1,0)

#---
#indicator of whether person is censored before their next visit
#this is used in the analysis when estimating the censoring weights
#C=1 if a person is censored before their next visit and 0 otherwise

dat.long=dat.long%>%group_by(id)%>%mutate(ones=1)%>%mutate(rownum=cumsum(ones))#row number within individual
dat.long=dat.long%>%group_by(id)%>%mutate(maxrow=max(rownum))#number of rows for each individual

dat.long$C=0
dat.long$C=ifelse(dat.long$maxrow==dat.long$rownum & dat.long$D.obs==0,1,dat.long$C)

