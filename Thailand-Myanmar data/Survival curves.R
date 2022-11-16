### Data analysis of multiple recurrence P. vivax data
# Survival curves


### Overview:
# 1. Load packages and data
# 2. Time to first recurrence: blood-stage vs primaquine + blood-stage treatment
# 3. Survival curves for time to 1st recurrence by drug and study
# 4. Survival curves for time to 2nd recurrence by drug and study
# 5. Survival curves for time to 2nd recurrence by time to first recurrence quartiles


#################################
### 1. Load packages and data ###
#################################

# load required packages:
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(splines)
library(R.matlab)

# load data:
load("./Data/Combined_Time_Event.RData")


######################################################################################
### 2. Time to first recurrence: blood-stage vs primaquine + blood-stage treatment ###
######################################################################################

### Survival curves for time to first recurrence: blood-stage treatment vs primaquine and blood-stage treatment (Fig. 1A)
tmin <- 0 # use tmin and tmax to plot only a part of the survival curve
tmax <- max(Combined_Time_Data$Time_to_event[which(Combined_Time_Data$episode==2)])+1

data_tmp <- Combined_Time_Data[which(Combined_Time_Data$episode==2 & Combined_Time_Data$Time_to_event>tmin),c(1,3,4,5)]
data_tmp$Censored <- 1-data_tmp$Censored
colnames(data_tmp) <- c("id","drug","time","event")
data_tmp$time <- data_tmp$time-tmin

# group patients by blood-stage treatment and primaquine & blood-stage treatment
data_tmp$drug[which(data_tmp$drug=="CHQ/PMQ")] <- "PMQ+"
data_tmp$drug[which(data_tmp$drug!="PMQ+")] <- "blood-stage"
fit_tmp <- survfit(Surv(time=data_tmp$time,event=data_tmp$event) ~ drug, data = data_tmp)
summary(fit_tmp)
ggsurvplot(fit_tmp, data = data_tmp, title="Time to first recurrence", pval = TRUE,break.time.by=100,
           # conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(), # for confidence intervals and risk table
           legend.title = "Treatment:", legend.labs = c("blood-stage treatment","primaquine"), legend="bottom", xlim=c(0,tmax-tmin))


### Incidence rate per patients at risk for blood-stage treatment vs primaquine and blood-stage treatment (Fig. 1C)
dt <- 7 # time steps (7 days for weekly incidence)
inc_bl <- cbind(seq(dt,max(data_tmp$time[data_tmp$event==1]),dt),matrix(0,length(seq(dt,max(data_tmp$time[data_tmp$event==1]),dt)))) # incidence for blood-stage treatment
inc_pmq <- cbind(seq(dt,max(data_tmp$time[data_tmp$event==1]),dt),matrix(0,length(seq(dt,max(data_tmp$time[data_tmp$event==1]),dt)))) # incidence for primaquine and blood-stage treatment
for(i in 1:dim(inc_bl)[1]){
  if(i==1){
    inc_bl[i,2] <- sum(data_tmp$time[data_tmp$drug=="blood-stage" & data_tmp$event==1]<=inc_bl[i,1])/sum(data_tmp$time[data_tmp$drug=="blood-stage"]>=0) # new cases/number of patients at risk
    inc_pmq[i,2] <- sum(data_tmp$time[data_tmp$drug=="PMQ+" & data_tmp$event==1]<=inc_pmq[i,1])/sum(data_tmp$time[data_tmp$drug=="PMQ+"]>=0) # new cases/number of patients at risk
  }else{
    inc_bl[i,2] <- sum(data_tmp$time[data_tmp$drug=="blood-stage" & data_tmp$event==1]>=inc_bl[i-1,1] & data_tmp$time[data_tmp$drug=="blood-stage" & data_tmp$event==1]<=inc_bl[i,1])/
      sum(data_tmp$time[data_tmp$drug=="blood-stage"]>=inc_bl[i-1,1]) # new cases/number of patients at risk
    inc_pmq[i,2] <- sum(data_tmp$time[data_tmp$drug=="PMQ+" & data_tmp$event==1]>=inc_pmq[i-1,1] & data_tmp$time[data_tmp$drug=="PMQ+" & data_tmp$event==1]<=inc_pmq[i,1])/
      sum(data_tmp$time[data_tmp$drug=="PMQ+"]>=inc_pmq[i-1,1]) # new cases/number of patients at risk
  }
}

# Plot points and fit a smooth curve through the points
# use the weekly incidence rate data as points
data_bl <- as.data.frame(inc_bl)
data_pmq <- as.data.frame(inc_pmq)
fit_bl <-lm(V2 ~ bs(V1,knots = seq(60,400,by=90)),data = data_bl)
# summary(fit_bl)
fit_pmq <-lm(V2 ~ bs(V1,knots = c(90,200,300,400)),data = data_pmq)
# summary(fit_pmq)

plot(inc_bl,pch=16,col="#F8766D",las=1,xlab="Time [days]",ylab="Weekly incidence rate",
     main="Weekly incidence rate per patients at risk",xlim=c(0,max(inc_bl[,1],inc_pmq[,1])))
points(1:400,predict(fit_bl,newdata = list(V1=1:400)),col="#F8766D",lwd=2,type="l")
points(inc_pmq,col="#00BFC4",pch=16)
points(1:400,predict(fit_pmq,newdata = list(V1=1:400)),col="#00BFC4",lwd=2,type="l")
legend("topright",legend=c("PMQ+ treatment", "blood-stage treatment"),col=c("#00BFC4","#F8766D"), lty = 1, lwd=2)


### Percent of blood-stage treated patients having a first infection for different time intervals
t_min <- 60
t_max <- 90
# number of patients with blood-stage treatment and no infection until t_min:
n_risk <- length(unique(Combined_Time_Data$patientid[which(Combined_Time_Data$episode==2 & 
                                                   Combined_Time_Data$arm_num%in%c("AS","CHQ") &
                                                   Combined_Time_Data$Time_to_event>t_min)]))
# number of patients with blood-stage treatment and first infection in (t_min,t_max]:
n_events <- length(unique(Combined_Time_Data$patientid[which(Combined_Time_Data$episode==2 & 
                                                     Combined_Time_Data$Censored==0 &
                                                     Combined_Time_Data$arm_num%in%c("AS","CHQ") &
                                                     Combined_Time_Data$Time_to_event>t_min &
                                                     Combined_Time_Data$Time_to_event<=t_max)]))
perc_event <- 100*n_events/n_risk


#######################################################################
### 3. Survival curves for time to 1st recurrence by drug and study ###
#######################################################################

data_tmp2 <- cbind(Combined_Time_Data,study=substr(Combined_Time_Data$patientid,1,3))
data_tmp2$arm_num[data_tmp2$PMQ_partner=="DP"] <- "DP/PMQ"
data_tmp2 <- cbind(data_tmp2,drug_study=paste(data_tmp2$arm_num,"_",data_tmp2$study,sep=""))
data_tmp2 <- data_tmp2[which(data_tmp2$episode==2),c(1,3,4,5,11,12)]
data_tmp2$Censored <- 1-data_tmp2$Censored
colnames(data_tmp2) <- c("id","drug","time","event","study","drug_study")
fit_tmp2 <- survfit(Surv(time=data_tmp2$time,event=data_tmp2$event) ~ drug_study, data = data_tmp2)
summary(fit_tmp2)

# Visualize without confidence region (Fig. S8):
ggsurvplot(fit_tmp2, data = data_tmp2, title="Time to first recurrence", pval = TRUE,
           legend.title = "Treatment:", legend.labs = c("AS","CHQ","CHQ/PMQ_BPD","CHQ/PMQ_VHX","DP/PMQ"), legend="bottom", xlim=c(0,400))

# Visualize with confidence region (confidence region used in Fig. S10, Fig. 3 and Fig. S11 first column):
ggsurvplot(fit_tmp2, data = data_tmp2, title="Time to first recurrence", pval = TRUE,
           conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=100,
           legend.title = "Treatment:", legend.labs = c("AS","CHQ","CHQ/PMQ_BPD","CHQ/PMQ_VHX","DP/PMQ"), legend="bottom", xlim=c(0,400))

# Survival curve data export to Matlab (for plotting together with model fits):
# surv_curve_data <- as.data.frame(cbind(summary(fit_tmp2)$time,summary(fit_tmp2)$surv,substr(summary(fit_tmp2)$strata,12,22)))
# names(surv_curve_data) <- c("time","surv","drug_study")
# surv_curve_data$drug_study <- factor(surv_curve_data$drug)
# surv_curve_data$time <- as.numeric(as.character(surv_curve_data$time))
# surv_curve_data$surv <- as.numeric(as.character(surv_curve_data$surv))
# writeMat("surv_curve_data_by_drug_and_study.mat", labpcexport = surv_curve_data)


#######################################################################
### 4. Survival curves for time to 2nd recurrence by drug and study ###
#######################################################################

data_tmp3 <- cbind(Combined_Time_Data,study=substr(Combined_Time_Data$patientid,1,3))
data_tmp3$arm_num[data_tmp3$PMQ_partner=="DP"] <- "DP/PMQ"
data_tmp3 <- cbind(data_tmp3,drug_study=paste(data_tmp3$arm_num,"_",data_tmp3$study,sep=""))
data_tmp3 <- data_tmp3[which(data_tmp3$episode==3),c(1,3,4,5,11,12)]
data_tmp3$Censored <- 1-data_tmp3$Censored
colnames(data_tmp3) <- c("id","drug","time","event","study","drug_study")
fit_tmp3 <- survfit(Surv(time=data_tmp3$time,event=data_tmp3$event) ~ drug_study, data = data_tmp3)
summary(fit_tmp3)

# Visualize (confidence region used in Fig. 3 and Fig. S11 second column):
ggsurvplot(fit_tmp3, data = data_tmp3, title="Time to second recurrence", pval = TRUE, 
           conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=100, # for confidence intervals and risk table
           legend.title = "Treatment:", legend.labs = c("AS","CHQ","CHQ/PMQ_BPD","CHQ/PMQ_VHX","DP/PMQ"), legend="bottom", xlim=c(0,400))


###########################################################################################
### 5. Survival curves for time to 2nd recurrence by time to first recurrence quartiles ###
###########################################################################################

# Fig. S7, Fig. 4C & D and Fig. S12 survival curves for the data:
data_tmp4 <- cbind(Combined_Time_Data,study=substr(Combined_Time_Data$patientid,1,3))
data_tmp4$arm_num[data_tmp4$PMQ_partner=="DP"] <- "DP/PMQ"
data_tmp4 <- cbind(data_tmp4,drug_study=paste(data_tmp4$arm_num,"_",data_tmp4$study,sep=""))
id <- data_tmp4$patientid[which(data_tmp4$episode==3)] # individuals with at least one recurrence
time1 <- data_tmp4$Time_to_event[data_tmp4$patientid%in%id & data_tmp4$episode==2] # time to first recurrence
event1 <- 1-data_tmp4$Censored[data_tmp4$patientid%in%id & data_tmp4$episode==2] # event at first recurrence
time2 <- data_tmp4$Time_to_event[data_tmp4$patientid%in%id & data_tmp4$episode==3] # time to second recurrence
event2 <- 1-data_tmp4$Censored[data_tmp4$patientid%in%id & data_tmp4$episode==3] # event at second recurrence
drug <- data_tmp4$arm_num[data_tmp4$patientid%in%id & data_tmp4$episode==2]
study <- data_tmp4$study[data_tmp4$patientid%in%id & data_tmp4$episode==2]
drug_study <- data_tmp4$drug_study[data_tmp4$patientid%in%id & data_tmp4$episode==2]
data_tmp4 <- data.frame(id=id,drug=drug,time1=as.numeric(time1),event1=as.numeric(event1),time2=as.numeric(time2),event2=as.numeric(event2),study,drug_study)

# group by time to first recurrence quartiles:
# for all data
data_tmp4_all <- as.data.frame(cbind(data_tmp4, time1_group = ifelse(data_tmp4$time1<summary(data_tmp4$time1)[2],paste("0-",summary(data_tmp4$time1)[2]-1,sep=""),
                                                                        ifelse(data_tmp4$time1<summary(data_tmp4$time1)[3],paste(summary(data_tmp4$time1)[2],"-",round(summary(data_tmp4$time1)[3]-1),sep=""),
                                                                               ifelse(data_tmp4$time1<round(summary(data_tmp4$time1)[5]),paste(summary(data_tmp4$time1)[3],"-",round(summary(data_tmp4$time1)[5]-1),sep=""),
                                                                                      paste(round(summary(data_tmp4$time1)[5]),"-",max(data_tmp4$time1),sep=""))))))
# for artesunate treatment (all in the VHX study)
data_tmp4_as <- data_tmp4[which(data_tmp4$study=="VHX" & data_tmp4$drug=="AS"),]
data_tmp4_as <- as.data.frame(cbind(data_tmp4_as, time1_group = ifelse(data_tmp4_as$time1<summary(data_tmp4_as$time1)[2],paste("0-",summary(data_tmp4_as$time1)[2]-1,sep=""),
                                                                     ifelse(data_tmp4_as$time1<summary(data_tmp4_as$time1)[3],paste(summary(data_tmp4_as$time1)[2],"-",round(summary(data_tmp4_as$time1)[3]-1),sep=""),
                                                                            ifelse(data_tmp4_as$time1<round(summary(data_tmp4_as$time1)[5]),paste(summary(data_tmp4_as$time1)[3],"-",round(summary(data_tmp4_as$time1)[5]-1),sep=""),
                                                                                   paste(round(summary(data_tmp4_as$time1)[5]),"-",max(data_tmp4_as$time1),sep=""))))))
# for chloroquine treatment (all in the VHX study)
data_tmp4_chq <- data_tmp4[which(data_tmp4$study=="VHX" & data_tmp4$drug=="CHQ"),]
data_tmp4_chq <- as.data.frame(cbind(data_tmp4_chq, time1_group = ifelse(data_tmp4_chq$time1<summary(data_tmp4_chq$time1)[2],paste("0-",summary(data_tmp4_chq$time1)[2]-1,sep=""),
                                                                       ifelse(data_tmp4_chq$time1<summary(data_tmp4_chq$time1)[3],paste(summary(data_tmp4_chq$time1)[2],"-",round(summary(data_tmp4_chq$time1)[3]-1),sep=""),
                                                                              ifelse(data_tmp4_chq$time1<round(summary(data_tmp4_chq$time1)[5]),paste(summary(data_tmp4_chq$time1)[3],"-",round(summary(data_tmp4_chq$time1)[5]-1),sep=""),
                                                                                     paste(round(summary(data_tmp4_chq$time1)[5]),"-",max(data_tmp4_chq$time1),sep=""))))))
# for chloroquine and primaquine treatment and the VHX study
data_tmp4_vhx_cpmq <- data_tmp4[which(data_tmp4$study=="VHX" & data_tmp4$drug=="CHQ/PMQ"),]
data_tmp4_vhx_cpmq <- as.data.frame(cbind(data_tmp4_vhx_cpmq, time1_group = ifelse(data_tmp4_vhx_cpmq$time1<summary(data_tmp4_vhx_cpmq$time1)[2],paste("0-",summary(data_tmp4_vhx_cpmq$time1)[2]-1,sep=""),
                                                                         ifelse(data_tmp4_vhx_cpmq$time1<summary(data_tmp4_vhx_cpmq$time1)[3],paste(summary(data_tmp4_vhx_cpmq$time1)[2],"-",round(summary(data_tmp4_vhx_cpmq$time1)[3]-1),sep=""),
                                                                                ifelse(data_tmp4_vhx_cpmq$time1<round(summary(data_tmp4_vhx_cpmq$time1)[5]),paste(summary(data_tmp4_vhx_cpmq$time1)[3],"-",round(summary(data_tmp4_vhx_cpmq$time1)[5])-2,sep=""),
                                                                                       paste(round(summary(data_tmp4_vhx_cpmq$time1)[5])-1,"-",max(data_tmp4_vhx_cpmq$time1),sep=""))))))
# for chloroquine and primaquine treatment and the BPD study
data_tmp4_bpd_cpmq <- data_tmp4[which(data_tmp4$study=="BPD" & data_tmp4$drug=="CHQ/PMQ"),]
data_tmp4_bpd_cpmq <- as.data.frame(cbind(data_tmp4_bpd_cpmq, time1_group = ifelse(data_tmp4_bpd_cpmq$time1<round(summary(data_tmp4_bpd_cpmq$time1)[2]),paste("0-",round(summary(data_tmp4_bpd_cpmq$time1)[2]),sep=""),
                                                                                   ifelse(data_tmp4_bpd_cpmq$time1<summary(data_tmp4_bpd_cpmq$time1)[3],paste(round(summary(data_tmp4_bpd_cpmq$time1)[2])+1,"-",round(summary(data_tmp4_bpd_cpmq$time1)[3]-1),sep=""),
                                                                                          ifelse(data_tmp4_bpd_cpmq$time1<round(summary(data_tmp4_bpd_cpmq$time1)[5]),paste(summary(data_tmp4_bpd_cpmq$time1)[3],"-",round(summary(data_tmp4_bpd_cpmq$time1)[5])-1,sep=""),
                                                                                                 paste(round(summary(data_tmp4_bpd_cpmq$time1)[5]),"-",max(data_tmp4_bpd_cpmq$time1),sep=""))))))
# for dihydroartemisinin-piperaquine and primaquine treatment (all in the BPD study)
data_tmp4_dpmq <- data_tmp4[which(data_tmp4$study=="BPD" & data_tmp4$drug=="DP/PMQ"),]
data_tmp4_dpmq <- as.data.frame(cbind(data_tmp4_dpmq, time1_group = ifelse(data_tmp4_dpmq$time1<summary(data_tmp4_dpmq$time1)[2],paste("0-",summary(data_tmp4_dpmq$time1)[2]-1,sep=""),
                                                                         ifelse(data_tmp4_dpmq$time1<round(summary(data_tmp4_dpmq$time1)[3]),paste(round(summary(data_tmp4_dpmq$time1)[2]),"-",round(summary(data_tmp4_dpmq$time1)[3])-1,sep=""),
                                                                                ifelse(data_tmp4_dpmq$time1<round(summary(data_tmp4_dpmq$time1)[5]),paste(round(summary(data_tmp4_dpmq$time1)[3]),"-",round(summary(data_tmp4_dpmq$time1)[5]-1),sep=""),
                                                                                       paste(round(summary(data_tmp4_dpmq$time1)[5])+1,"-",max(data_tmp4_dpmq$time1),sep=""))))))

# Plot the survival curves:
# for all data
fit_tmp4_all <- survfit(Surv(time=data_tmp4_all$time2,event=data_tmp4_all$event2) ~ data_tmp4_all$time1_group, data = data_tmp4_all)
summary(fit_tmp4_all)
ggsurvplot(fit_tmp4_all, data = data_tmp4_all, title="Time from first to second recurrence (all data)", pval = TRUE, 
           conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=100, # for confidence intervals and risk table
           legend.title = "Time to first recurrence:", legend="bottom")
# for artesunate treatment (all in the VHX study)
fit_tmp4_as <- survfit(Surv(time=data_tmp4_as$time2,event=data_tmp4_as$event2) ~ data_tmp4_as$time1_group, data = data_tmp4_as)
summary(fit_tmp4_as)
ggsurvplot(fit_tmp4_as, data = data_tmp4_as, title="Time from first to second recurrence (AS treatment)", pval = TRUE, 
           conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=100, # for confidence intervals and risk table
           legend.title = "Time to first recurrence:", legend="bottom")
# for chloroquine treatment (all in the VHX study)
fit_tmp4_chq <- survfit(Surv(time=data_tmp4_chq$time2,event=data_tmp4_chq$event2) ~ data_tmp4_chq$time1_group, data = data_tmp4_chq)
summary(fit_tmp4_chq)
ggsurvplot(fit_tmp4_chq, data = data_tmp4_chq, title="Time from first to second recurrence (CHQ treatment)", pval = TRUE, 
           conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=100, # for confidence intervals and risk table
           legend.title = "Time to first recurrence:", legend="bottom")
# for chloroquine and primaquine treatment and the VHX study
fit_tmp4_vhx_cpmq <- survfit(Surv(time=data_tmp4_vhx_cpmq$time2,event=data_tmp4_vhx_cpmq$event2) ~ data_tmp4_vhx_cpmq$time1_group, data = data_tmp4_vhx_cpmq)
summary(fit_tmp4_vhx_cpmq)
ggsurvplot(fit_tmp4_vhx_cpmq, data = data_tmp4_vhx_cpmq, title="Time from first to second recurrence (CHQ/PMQ treatment, VHX study)", pval = TRUE, 
           conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=100, # for confidence intervals and risk table
           legend.title = "Time to first recurrence:", legend="bottom")
# for chloroquine and primaquine treatment and the BPD study
fit_tmp4_bpd_cpmq <- survfit(Surv(time=data_tmp4_bpd_cpmq$time2,event=data_tmp4_bpd_cpmq$event2) ~ data_tmp4_bpd_cpmq$time1_group, data = data_tmp4_bpd_cpmq)
summary(fit_tmp4_bpd_cpmq)
ggsurvplot(fit_tmp4_bpd_cpmq, data = data_tmp4_bpd_cpmq, title="Time from first to second recurrence (CHQ/PMQ treatment, BPD study)", pval = TRUE, 
           conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=100, # for confidence intervals and risk table
           legend.title = "Time to first recurrence:", legend="bottom")
# for dihydroartemisinin-piperaquine and primaquine treatment (all in the BPD study)
fit_tmp4_dpmq <- survfit(Surv(time=data_tmp4_dpmq$time2,event=data_tmp4_dpmq$event2) ~ data_tmp4_dpmq$time1_group, data = data_tmp4_dpmq)
summary(fit_tmp4_dpmq)
ggsurvplot(fit_tmp4_dpmq, data = data_tmp4_dpmq, title="Time from first to second recurrence (DP/PMQ treatment, BPD study)", pval = TRUE, 
           conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=100, # for confidence intervals and risk table
           legend.title = "Time to first recurrence:", legend="bottom")


