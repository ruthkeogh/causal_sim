#------------------------------------------------------
#------------------------------------------------------
# dat_sim_TRUE_VALUES_1.R
# - Generates data from a large trial to obtain (estimated) true values of 
# the cumulative coefficients of the correctly-specified MSM, 
# and of the causal estimands of interest (survival probabilities under the treatment regimes 'always treated' and 'never treated'). 
# See Section 6.2. of the paper
# Created by Ruth Keogh
#------------------------------------------------------
#------------------------------------------------------

#----
#number of individuals per treatment pattern 

n=1000

#----
#number of visits (K+1)

n.visit=5

#----
#treatment patterns
#there are 32=2^5 possible treatment patterns across 5 visits

A.patterns=matrix(nrow=32,ncol=5)
A.patterns[1,]=c(1,1,1,1,1)
A.patterns[2,]=c(0,0,0,0,0)

A.patterns[3,]=c(1,0,0,0,0)
A.patterns[4,]=c(0,1,0,0,0)
A.patterns[5,]=c(0,0,1,0,0)
A.patterns[6,]=c(0,0,0,1,0)
A.patterns[7,]=c(0,0,0,0,1)

A.patterns[8,]=c(1,1,0,0,0)
A.patterns[9,]=c(1,0,1,0,0)
A.patterns[10,]=c(1,0,0,1,0)
A.patterns[11,]=c(1,0,0,0,1)
A.patterns[12,]=c(0,1,1,0,0)
A.patterns[13,]=c(0,1,0,1,0)
A.patterns[14,]=c(0,1,0,0,1)
A.patterns[15,]=c(0,0,1,1,0)
A.patterns[16,]=c(0,0,1,0,1)
A.patterns[17,]=c(0,0,0,1,1)

A.patterns[18,]=c(1,1,1,0,0)
A.patterns[19,]=c(1,1,0,1,0)
A.patterns[20,]=c(1,1,0,0,1)
A.patterns[21,]=c(1,0,1,1,0)
A.patterns[22,]=c(1,0,1,0,1)
A.patterns[23,]=c(1,0,0,1,1)
A.patterns[24,]=c(0,1,1,1,0)
A.patterns[25,]=c(0,1,1,0,1)
A.patterns[26,]=c(0,1,0,1,1)
A.patterns[27,]=c(0,0,1,1,1)

A.patterns[28,]=c(0,1,1,1,1)
A.patterns[29,]=c(1,0,1,1,1)
A.patterns[30,]=c(1,1,0,1,1)
A.patterns[31,]=c(1,1,1,0,1)
A.patterns[32,]=c(1,1,1,1,0)

#-------------------
#data generation: a data set is created under each treatment pattern, and then they are combined
#-------------------

#----
#make L0 the same in all treatment pattern groups
U=rnorm(n,0,0.1)
L0.common=rnorm(n,0+U,1)

#----
#for each treatment pattern

for(pattern in 1:32){
  
  #generate L at each visit
  L=matrix(nrow=n,ncol=n.visit)
  A=matrix(nrow=n,ncol=n.visit)
  
  L[,1]=L0.common
  A[,1]=A.patterns[pattern,1]
  for(k in 2:n.visit){
    A[,k]=A.patterns[pattern,k]
    L[,k]=rnorm(n,0.8*L[,k-1]-A[,k-1]+0.1*(k-1)+U,1)
  }
  
  #generate event times T.obs, and event indicators D.obs
  T.obs=rep(NA,n)
  
  for(v in 1:5){
    u.t=runif(n,0,1)
    haz=0.7-0.2*A[,v]+0.05*L[,v]+0.05*U
    new.t=-log(u.t)/haz
    T.obs=ifelse(is.na(T.obs) & new.t<1 & haz>0,v-1+new.t,T.obs) #the haz>0 is just used to deal with tiny possibility (under this data generating mechanism) the hazard could go negative. 
  }
  D.obs=ifelse(is.na(T.obs),0,1)
  T.obs=ifelse(is.na(T.obs),5,T.obs)
  
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
  
  dat=data.frame(id=1:n,T.obs,D.obs,L.dat,A.dat,Alag1.dat,Alag2.dat,Alag3.dat,Alag4.dat)
  
  dat.long=reshape(data = dat,varying=c(paste0("A.",0:4),paste0("Alag1.",0:4),paste0("Alag2.",0:4),paste0("Alag3.",0:4),paste0("Alag4.",0:4),paste0("L.",0:4)),direction="long",idvar="id")
  
  dat.long=dat.long[order(dat.long$id,dat.long$time),]
  
  dat.long$time.stop=dat.long$time+1
  
  dat.long=dat.long[dat.long$time<dat.long$T.obs,]
  
  dat.long$time.stop=ifelse(dat.long$time.stop>dat.long$T.obs,dat.long$T.obs,dat.long$time.stop)
  
  dat.long$event=ifelse(dat.long$time.stop==dat.long$T.obs & dat.long$D.obs==1,1,0)
  
  if(pattern==1){
    dat.long.combine=dat.long
  }else if(pattern>1){
    dat.long.combine=rbind(dat.long.combine,dat.long)
  }
}

#-------------------
#data set used to obtain true values of model parameters and causal estimands
#-------------------

dat.long=dat.long.combine


