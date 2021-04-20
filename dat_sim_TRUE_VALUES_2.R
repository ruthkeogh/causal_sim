#------------------------------------------------------
#------------------------------------------------------
# dat_sim_TRUE_VALUES_2.R
# - Generates data from a large trial focusing on only two treatment regimes: 'always treated' and 'never treated'.
# - Generates data from a large trial to obtain (estimated) true values of
# the causal estimands of interest (survival probabilities under the treatment regimes 'always treated' and 'never treated'). 
# - Unlike dat_sim_TRUE_VALUES_1.R, this file only simulates data under two treatment regimes: 'always treated' and 'never treated'.
# See Section 6.2. of the paper
# Created by Ruth Keogh
#------------------------------------------------------
#------------------------------------------------------

#----
#number of individuals per treatment pattern (focusing only on two patterns: 'always treated' and 'never treated')

n=10000

#----
#number of visits (K+1)

n.visit=5


#----
#treatment patterns: 'always treated' and 'never treated'

A.patterns=matrix(nrow=2,ncol=5)
A.patterns[1,]=c(1,1,1,1,1)
A.patterns[2,]=c(0,0,0,0,0)

#-------------------
#data generation: a data set is created under each treatment pattern
#-------------------

#----
#make L0 the same in both treatment pattern groups
U=rnorm(n,0,0.1)
L0.common=rnorm(n,0+U,1)

#----
#for each treatment pattern
for(pattern in 1:2){
  
  L=matrix(nrow=n,ncol=n.visit)
  A=matrix(nrow=n,ncol=n.visit)
  
  #generate L at each visit
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
  
  if(pattern==1){
    dat.1=data.frame(id=1:n,A=1,D.obs,T.obs) #data set for the 'always treated'
  }
  if(pattern==2){
    dat.0=data.frame(id=(n+1):(2*n),A=0,D.obs,T.obs) #data set for the 'never treated'
  }
  
}


