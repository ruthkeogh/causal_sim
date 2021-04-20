#------------------------------------------------------
#------------------------------------------------------
# sim_results.R
# - Summarises the simulation results to create Tables 1 and 2 and Figures 2 and 3
# Created by Ruth Keogh
#------------------------------------------------------
#------------------------------------------------------

#------------------
#time increments: 
#this is a grid of times from 0 to 5 at which the cumulative coefficients and survival probabilities are evaluated
#------------------

t.hor=seq(0,5,0.01)

pos.t1=which(t.hor==1)
pos.t2=which(t.hor==2)
pos.t3=which(t.hor==3)
pos.t4=which(t.hor==4)
pos.t5=which(t.hor==5)

#------------------
#cumulative coefficients at visits k=1,2,3,4,5
#See Table 1 in the paper
#------------------

#----
#true values

est.coef.0.true=coef_0_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
est.coef.A.true=coef_A_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
est.coef.Alag1.true=coef_Alag1_true[c(pos.t2,pos.t3,pos.t4,pos.t5)]
est.coef.Alag2.true=coef_Alag2_true[c(pos.t3,pos.t4,pos.t5)]
est.coef.Alag3.true=coef_Alag3_true[c(pos.t4,pos.t5)]
est.coef.Alag4.true=coef_Alag4_true[pos.t5]

#----
#mean of estimates using MSM-IPTW

est.coef.0.msm=colMeans(coef_0)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
est.coef.A.msm=colMeans(coef_A)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
est.coef.Alag1.msm=colMeans(coef_Alag1)[c(pos.t2,pos.t3,pos.t4,pos.t5)]
est.coef.Alag2.msm=colMeans(coef_Alag2)[c(pos.t3,pos.t4,pos.t5)]
est.coef.Alag3.msm=colMeans(coef_Alag3)[c(pos.t4,pos.t5)]
est.coef.Alag4.msm=colMeans(coef_Alag4)[c(pos.t5)]

#----
#SD of estimates (empirical standard errors) using MSM-IPTW

empsd.coef.0.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_0[,x])})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
empsd.coef.A.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_A[,x])})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
empsd.coef.Alag1.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag1[,x])})[c(pos.t2,pos.t3,pos.t4,pos.t5)]
empsd.coef.Alag2.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag2[,x])})[c(pos.t3,pos.t4,pos.t5)]
empsd.coef.Alag3.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag3[,x])})[c(pos.t4,pos.t5)]
empsd.coef.Alag4.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag4[,x])})[c(pos.t5)]

#----
#bias in estimates obtained using MSM-IPTW

bias.coef.0.msm=colMeans(coef_0)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]-coef_0_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
bias.coef.A.msm=colMeans(coef_A)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]-coef_A_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
bias.coef.Alag1.msm=colMeans(coef_Alag1)[c(pos.t2,pos.t3,pos.t4,pos.t5)]-coef_Alag1_true[c(pos.t2,pos.t3,pos.t4,pos.t5)]
bias.coef.Alag2.msm=colMeans(coef_Alag2)[c(pos.t3,pos.t4,pos.t5)]-coef_Alag2_true[c(pos.t3,pos.t4,pos.t5)]
bias.coef.Alag3.msm=colMeans(coef_Alag3)[c(pos.t4,pos.t5)]-coef_Alag3_true[c(pos.t4,pos.t5)]
bias.coef.Alag4.msm=colMeans(coef_Alag4)[c(pos.t5)]-coef_Alag4_true[c(pos.t5)]

#----
#Monte Carlo standard errors for estimates obtained using MSM-IPTW

mcbias.coef.0.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_0[,x]-coef_0_true[x])/sqrt(n.sim)})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
mcbias.coef.A.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_A[,x]-coef_A_true[x])/sqrt(n.sim)})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
mcbias.coef.Alag1.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag1[,x]-coef_Alag1_true[x])/sqrt(n.sim)})[c(pos.t2,pos.t3,pos.t4,pos.t5)]
mcbias.coef.Alag2.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag2[,x]-coef_Alag2_true[x])/sqrt(n.sim)})[c(pos.t3,pos.t4,pos.t5)]
mcbias.coef.Alag3.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag3[,x]-coef_Alag3_true[x])/sqrt(n.sim)})[c(pos.t4,pos.t5)]
mcbias.coef.Alag4.msm=sapply(1:length(t.hor),FUN=function(x){sd(coef_Alag4[,x]-coef_Alag4_true[x])/sqrt(n.sim)})[c(pos.t5)]

#----
#table of results: as in Table 1

#construct the table
table.coef.0=data.frame(Cumulative_coefficient="alpha_0(t)",Time=1:5,
                        True_value=dp3(est.coef.0.true),
                        Mean_estimate__Empirical_SE=table.func(est.coef.0.msm,empsd.coef.0.msm),
                        Bias__MonteCarlo_SE=table.func(bias.coef.0.msm,mcbias.coef.0.msm))

table.coef.A=data.frame(Cumulative_coefficient="alpha_A,0(t)",Time=1:5,
                        True_value=dp3(est.coef.A.true),
                        Mean_estimate__Empirical_SE=table.func(est.coef.A.msm,empsd.coef.A.msm),
                        Bias__MonteCarlo_SE=table.func(bias.coef.A.msm,mcbias.coef.A.msm))

table.coef.Alag1=data.frame(Cumulative_coefficient="alpha_A,1(t)",Time=2:5,
                            True_value=dp3(est.coef.Alag1.true),
                            Mean_estimate__Empirical_SE=table.func(est.coef.Alag1.msm,empsd.coef.Alag1.msm)[1:4],
                            Bias__MonteCarlo_SE=table.func(bias.coef.Alag1.msm,mcbias.coef.Alag1.msm)[1:4])

table.coef.Alag2=data.frame(Cumulative_coefficient="alpha_A,2(t)",Time=3:5,
                            True_value=dp3(est.coef.Alag2.true),
                            Mean_estimate__Empirical_SE=table.func(est.coef.Alag2.msm,empsd.coef.Alag2.msm)[1:3],
                            Bias__MonteCarlo_SE=table.func(bias.coef.Alag2.msm,mcbias.coef.Alag2.msm)[1:3])

table.coef.Alag3=data.frame(Cumulative_coefficient="alpha_A,3(t)",Time=4:5,
                            True_value=dp3(est.coef.Alag3.true),
                            Mean_estimate__Empirical_SE=table.func(est.coef.Alag3.msm,empsd.coef.Alag3.msm)[1:2],
                            Bias__MonteCarlo_SE=table.func(bias.coef.Alag3.msm,mcbias.coef.Alag3.msm)[1:2])

table.coef.Alag4=data.frame(Cumulative_coefficient="alpha_A,4(t)",Time=5,
                            True_value=dp3(est.coef.Alag4.true),
                            Mean_estimate__Empirical_SE=table.func(est.coef.Alag4.msm,empsd.coef.Alag4.msm)[1],
                            Bias__MonteCarlo_SE=table.func(bias.coef.Alag4.msm,mcbias.coef.Alag4.msm)[1])

table1=rbind(table.coef.0,table.coef.A,table.coef.Alag1,table.coef.Alag2,table.coef.Alag3,table.coef.Alag4)

#save table as a csv file

write.table(table1,"results/table1.csv",sep=",",row.names = F)

#------------------
#survival probabilities at visits k=1,2,3,4,5
#See Table 2 in the paper
#------------------

#----
#true values

est.surv1.true=surv1_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
est.surv0.true=surv0_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]

#----
#mean of estimates using MSM-IPTW

est.surv1.msm=colMeans(surv1)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
est.surv0.msm=colMeans(surv0)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]

#----
#SD of estimates (empirical standard errors) using MSM-IPTW

empsd.surv1.msm=sapply(1:length(t.hor),FUN=function(x){sd(surv1[,x])})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
empsd.surv0.msm=sapply(1:length(t.hor),FUN=function(x){sd(surv0[,x])})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]

#----
#bias in estimates obtained using MSM-IPTW

bias.surv1.msm=colMeans(surv1)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]-surv1_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
bias.surv0.msm=colMeans(surv0)[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]-surv0_true[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]

#----
#MC standard errors for estimates obtained using MSM-IPTW

mcbias.surv1.msm=sapply(1:length(t.hor),FUN=function(x){sd(surv1[,x]-surv1_true[x])/sqrt(n.sim)})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]
mcbias.surv0.msm=sapply(1:length(t.hor),FUN=function(x){sd(surv0[,x]-surv0_true[x])/sqrt(n.sim)})[c(pos.t1,pos.t2,pos.t3,pos.t4,pos.t5)]

#----
#latex tables for results: as in Table 2

table.surv.never.treated=data.frame(Treatment_regime="Never treated",Time=1:5,
                                    True_value=dp3(est.surv0.true),
                                    Mean_estimate__Empirical_SE=table.func(est.surv0.msm,empsd.surv0.msm),
                                    Bias__MonteCarlo_SE=table.func(bias.surv0.msm,mcbias.surv0.msm))

table.surv.always.treated=data.frame(Treatment_regime="Always treated",Time=1:5,
                                     True_value=dp3(est.surv1.true),
                                     Mean_estimate__Empirical_SE=table.func(est.surv1.msm,empsd.surv1.msm),
                                     Bias__MonteCarlo_SE=table.func(bias.surv1.msm,mcbias.surv1.msm))

table2=rbind(table.surv.never.treated,table.surv.always.treated)

#save table as a csv file

write.table(table2,"results/table2.csv",sep=",",row.names = F)

#------------------
#plot of cumulative coefficients
#See Figure 2 in the paper
#------------------

#colours used in the plot
cols=c("true"="#000000","msmiptw"="#999999")

#------
#create data frame containing true values and mean of MSM-IPTW estimates

dat.coef=data.frame(time=t.hor,coef_0=colMeans(coef_0),coef_0_true,
                    coef_A=colMeans(coef_A),coef_A_true,
                    coef_Alag1=colMeans(coef_Alag1),coef_Alag1_true,
                    coef_Alag2=colMeans(coef_Alag2),coef_Alag2_true,
                    coef_Alag3=colMeans(coef_Alag3),coef_Alag3_true,
                    coef_Alag4=colMeans(coef_Alag4),coef_Alag4_true)

#---------
#alpha_0

#create data frame containing MSM-IPTW estimates from each simulated data set
#this is used to add the faint lines on the plot
dat.coef0.lines=as.matrix(coef_0)
dat.coef0.lines=as.vector(t(dat.coef0.lines))
dat.coef0.lines=data.frame(dat.coef0.lines)
names(dat.coef0.lines)="coef"
dat.coef0.lines$time=rep(t.hor,n.sim)
dat.coef0.lines$sim=rep(1:n.sim,each=length(t.hor))

#create the plot
p0 <- ggplot(data = dat.coef0.lines, aes(x = time, y = coef, group = sim),alpha=0.01)
p0=p0+geom_line(size=0.0001,colour="gray94")
p0=p0+addlinetoplot(dat.coef,"time","coef_0_true",vcol='"true"')+
  addlinetoplot(dat.coef,"time","coef_0",vcol='"msmiptw"')+
  ylab(expression(paste("Cumulative coefficient ","  ",tilde(alpha)[0](t))))+xlab("Time")+
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(0,4,1),limits=c(0,4))+
  scale_colour_manual(NULL,values=cols,labels=c(true="True",msmiptw="MSM-IPTW"),breaks=c("true","msmiptw"))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),legend.text=element_text(size = 12))+
  theme(legend.position=c(0.35,0.8))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.key.width = unit(2, "cm"))+
  theme(legend.key = element_rect(fill = NA))

#---------
#alpha_A

#create data frame containing MSM-IPTW estimates from each simulated data set
#this is used to add the faint lines on the plot
dat.coefA.lines=as.matrix(coef_A)
dat.coefA.lines=as.vector(t(dat.coefA.lines))
dat.coefA.lines=data.frame(dat.coefA.lines)
names(dat.coefA.lines)="coef"
dat.coefA.lines$time=rep(t.hor,n.sim)
dat.coefA.lines$sim=rep(1:n.sim,each=length(t.hor))

#create the plot
pA <- ggplot(data = dat.coefA.lines, aes(x = time, y = coef, group = sim),alpha=0.01)
pA=pA+geom_line(size=0.0001,colour="gray94")
pA=pA+addlinetoplot(dat.coef,"time","coef_A_true",vcol='"true"')+
  addlinetoplot(dat.coef,"time","coef_A",vcol='"msmiptw"')+
  ylab(expression(paste("Cumulative coefficient ","  ",tilde(alpha)["A0"](t))))+xlab("Time")+
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(-1.6,0,0.2),limits=c(-1.6,0.1))+
  scale_colour_manual(NULL,values=cols,labels=c(true="True",msmiptw="MSM-IPTW"),breaks=c("true","msmiptw"))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),legend.text=element_text(size = 12))+
  theme(legend.position=c(0.35,0.2))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.key.width = unit(2, "cm"))+
  theme(legend.key = element_rect(fill = NA))

#---------
#alpha_Alag1

#create data frame containing MSM-IPTW estimates from each simulated data set
#this is used to add the faint lines on the plot
dat.coefAlag1.lines=as.matrix(coef_Alag1)
dat.coefAlag1.lines=as.vector(t(dat.coefAlag1.lines))
dat.coefAlag1.lines=data.frame(dat.coefAlag1.lines)
names(dat.coefAlag1.lines)="coef"
dat.coefAlag1.lines$time=rep(t.hor,n.sim)
dat.coefAlag1.lines$sim=rep(1:n.sim,each=length(t.hor))

#create the plot
pAlag1<- ggplot(data = dat.coefAlag1.lines, aes(x = time, y = coef, group = sim),alpha=0.01)
pAlag1=pAlag1+geom_line(size=0.0001,colour="gray94")
pAlag1=pAlag1+addlinetoplot(dat.coef,"time","coef_Alag1_true",vcol='"true"')+
  addlinetoplot(dat.coef,"time","coef_Alag1",vcol='"msmiptw"')+
  ylab(expression(paste("Cumulative coefficient ","  ",tilde(alpha)["A1"](t))))+xlab("Time")+
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(-1,1,0.2),limits=c(-1,1))+
  scale_colour_manual(NULL,values=cols,labels=c(true="True",msmiptw="MSM-IPTW"),breaks=c("true","msmiptw"))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),legend.text=element_text(size = 12))+
  theme(legend.position=c(0.35,0.2))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.key.width = unit(2, "cm"))+
  theme(legend.key = element_rect(fill = NA))

#---------
#alpha_Alag2

#create data frame containing MSM-IPTW estimates from each simulated data set
#this is used to add the faint lines on the plot
dat.coefAlag2.lines=as.matrix(coef_Alag2)
dat.coefAlag2.lines=as.vector(t(dat.coefAlag2.lines))
dat.coefAlag2.lines=data.frame(dat.coefAlag2.lines)
names(dat.coefAlag2.lines)="coef"
dat.coefAlag2.lines$time=rep(t.hor,n.sim)
dat.coefAlag2.lines$sim=rep(1:n.sim,each=length(t.hor))

#create the plot
pAlag2 <- ggplot(data = dat.coefAlag2.lines, aes(x = time, y = coef, group = sim),alpha=0.01)
pAlag2=pAlag2+geom_line(size=0.0001,colour="gray94")
pAlag2=pAlag2+addlinetoplot(dat.coef,"time","coef_Alag2_true",vcol='"true"')+
  addlinetoplot(dat.coef,"time","coef_Alag2",vcol='"msmiptw"')+
  ylab(expression(paste("Cumulative coefficient ","  ",tilde(alpha)["A2"](t))))+xlab("Time")+
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(-1,1,0.2),limits=c(-1,1))+
  scale_colour_manual(NULL,values=cols,labels=c(true="True",msmiptw="MSM-IPTW"),breaks=c("true","msmiptw"))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),legend.text=element_text(size = 12))+
  theme(legend.position=c(0.35,0.2))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.key.width = unit(2, "cm"))+
  theme(legend.key = element_rect(fill = NA))

#---------
#alpha_Alag3

#create data frame containing MSM-IPTW estimates from each simulated data set
#this is used to add the faint lines on the plot
dat.coefAlag3.lines=as.matrix(coef_Alag3)
dat.coefAlag3.lines=as.vector(t(dat.coefAlag3.lines))
dat.coefAlag3.lines=data.frame(dat.coefAlag3.lines)
names(dat.coefAlag3.lines)="coef"
dat.coefAlag3.lines$time=rep(t.hor,n.sim)
dat.coefAlag3.lines$sim=rep(1:n.sim,each=length(t.hor))

#create the plot
pAlag3 <- ggplot(data = dat.coefAlag3.lines, aes(x = time, y = coef, group = sim),alpha=0.01)
pAlag3=pAlag3+geom_line(size=0.0001,colour="gray94")
pAlag3=pAlag3+addlinetoplot(dat.coef,"time","coef_Alag3_true",vcol='"true"')+
  addlinetoplot(dat.coef,"time","coef_Alag3",vcol='"msmiptw"')+
  ylab(expression(paste("Cumulative coefficient ","  ",tilde(alpha)["A3"](t))))+xlab("Time")+
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(-1,1,0.2),limits=c(-1,1))+
  scale_colour_manual(NULL,values=cols,labels=c(true="True",msmiptw="MSM-IPTW"),breaks=c("true","msmiptw"))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),legend.text=element_text(size = 12))+
  theme(legend.position=c(0.35,0.2))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.key.width = unit(2, "cm"))+
  theme(legend.key = element_rect(fill = NA))

#---------
#alpha_Alag4

#create data frame containing MSM-IPTW estimates from each simulated data set
#this is used to add the faint lines on the plot
dat.coefAlag4.lines=as.matrix(coef_Alag4)
dat.coefAlag4.lines=as.vector(t(dat.coefAlag4.lines))
dat.coefAlag4.lines=data.frame(dat.coefAlag4.lines)
names(dat.coefAlag4.lines)="coef"
dat.coefAlag4.lines$time=rep(t.hor,n.sim)
dat.coefAlag4.lines$sim=rep(1:n.sim,each=length(t.hor))

#create the plot
pAlag4 <- ggplot(data = dat.coefAlag4.lines, aes(x = time, y = coef, group = sim),alpha=0.01)
pAlag4=pAlag4+geom_line(size=0.0001,colour="gray94")
pAlag4=pAlag4+addlinetoplot(dat.coef,"time","coef_Alag4_true",vcol='"true"')+
  addlinetoplot(dat.coef,"time","coef_Alag4",vcol='"msmiptw"')+
  ylab(expression(paste("Cumulative coefficient ","  ",tilde(alpha)["A4"](t))))+xlab("Time")+
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(-1,1,0.2),limits=c(-1,1))+
  scale_colour_manual(NULL,values=cols,labels=c(true="True",msmiptw="MSM-IPTW"),breaks=c("true","msmiptw"))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),legend.text=element_text(size = 12))+
  theme(legend.position=c(0.35,0.2))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.key.width = unit(2, "cm"))+
  theme(legend.key = element_rect(fill = NA))

#save subplots in figure 2 as png files
ggsave(p0,filename="results/figure2_1.png") #row 1, left in figure 2
ggsave(pA,filename="results/figure2_2.png") #row 1, right in figure 2
ggsave(pAlag1,filename="results/figure2_3.png") #row 2, left in figure 2
ggsave(pAlag2,filename="results/figure2_4.png") #row 2, right in figure 2
ggsave(pAlag3,filename="results/figure2_5.png") #row 3, left in figure 2
ggsave(pAlag4,filename="results/figure2_6.png") #row 3, right in figure 2

#the files can be combined using grid.arrange (for example),but this will look different depending on the size of your screen
#grid.arrange(p0,pA,pAlag1,pAlag2,pAlag3,pAlag4,ncol=2)

#------------------
#plot of survival curves
#See Figure 3 in the paper
#------------------

#colours used in the plot
cols.2=c("true0"="#000000","true1"="#000000","msmiptw0"="#999999","msmiptw1"="#999999")

#linetypes used in the plot
lines=c("true0"="dashed","true1"="solid","msmiptw0"="dashed","msmiptw1"="solid")

#------
#create data frame containing true values and mean of MSM-IPTW estimates

dat.surv=data.frame(time=t.hor,surv1=colMeans(surv1),surv0=colMeans(surv0),surv1_true,surv0_true)

#create data frames containing MSM-IPTW estimates from each simulated data set
#these are used to add the faint lines on the plot
dat.surv1.lines=as.matrix(surv1)
dat.surv1.lines=as.vector(t(dat.surv1.lines))
dat.surv1.lines=data.frame(dat.surv1.lines)
names(dat.surv1.lines)="surv"
dat.surv1.lines$time=rep(t.hor,n.sim)
dat.surv1.lines$sim=rep(1:n.sim,each=length(t.hor))

dat.surv0.lines=as.matrix(surv0)
dat.surv0.lines=as.vector(t(dat.surv0.lines))
dat.surv0.lines=data.frame(dat.surv0.lines)
names(dat.surv0.lines)="surv"
dat.surv0.lines$time=rep(t.hor,n.sim)
dat.surv0.lines$sim=rep((n.sim+1):(2*n.sim),each=length(t.hor))

dat.lines=data.frame(rbind(dat.surv1.lines,dat.surv0.lines))

#create the plot
psurv <- ggplot(data = dat.lines, aes(x = time, y = surv, group = sim),alpha=0.01)
psurv=psurv+geom_line(size=0.0001,colour="gray94")
psurv=psurv+addlinetoplot.2(dat.surv,"time","surv0_true",vcol='"true0"',vline='"true0"')+
  addlinetoplot.2(dat.surv,"time","surv1_true",vcol='"true1"',vline='"true1"')+
  addlinetoplot.2(dat.surv,"time","surv0",vcol='"msmiptw0"',vline='"msmiptw0"')+
  addlinetoplot.2(dat.surv,"time","surv1",vcol='"msmiptw1"',vline='"msmiptw1"')+
  ylab("Survival probability")+xlab("Time")+
  scale_x_continuous(breaks=seq(0,5,1),limits=c(0,5))+
  scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1))+
  scale_colour_manual(NULL,values=cols.2,labels=c(true0="True: Never treated",true1="True: Always treated",msmiptw0="MSM-IPTW: Never treated",msmiptw1="MSM-IPTW: Always treated"),breaks=c("true0","true1","msmiptw0","msmiptw1"))+
  scale_linetype_manual(NULL,values=lines,labels=c(true0="True: Never treated",true1="True: Always treated",msmiptw0="MSM-IPTW: Never treated",msmiptw1="MSM-IPTW: Always treated"),breaks=c("true0","true1","msmiptw0","msmiptw1"))+
  theme(axis.text.x = element_text(size = 12),axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12),legend.text=element_text(size = 12))+
  theme(legend.position=c(0.6,0.7))+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.background = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.key.width = unit(2, "cm"))+
  theme(legend.key = element_rect(fill = NA))

#save figure 3 as png files
ggsave(psurv,filename="results/figure3.png")
