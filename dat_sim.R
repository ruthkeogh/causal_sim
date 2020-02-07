#------------------------------------------------------
#------------------------------------------------------
# dat_sim.R
# - Simulates longitudinal and time-to-event data using the algorithm described in Section 5. 
#------------------------------------------------------
#------------------------------------------------------

#----
#sample size 

n=5000

#----
#number of visits (K+1)

n.visit=5

#----
#expit function

expit=function(x){exp(x)/(1+exp(x))}

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

dat=data.frame(id=1:n,T.obs,D.obs,A.dat,Alag1.dat,Alag2.dat,Alag3.dat,Alag4.dat,L.dat,U)

dat.long=reshape(data = dat,varying=c(paste0("A.",0:4),paste0("Alag1.",0:4),paste0("Alag2.",0:4),paste0("Alag3.",0:4),paste0("Alag4.",0:4),paste0("L.",0:4)),direction="long",idvar="id")

dat.long=dat.long[order(dat.long$id,dat.long$time),]

dat.long$time.stop=dat.long$time+1

dat.long=dat.long[dat.long$time<dat.long$T.obs,]

dat.long$time.stop=ifelse(dat.long$time.stop>dat.long$T.obs,dat.long$T.obs,dat.long$time.stop)

dat.long$event=ifelse(dat.long$time.stop==dat.long$T.obs & dat.long$D.obs==1,1,0)
