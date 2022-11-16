### Data analysis of P. vivax data from Robinson et al. (2015)


### Overview:
# 1. Load packages and data, combine data from different tables
# 2. Survival curves
# 3. Weekly incidence rate
# 4. Transmission and relapse rate (from model fits)


#####################################################################
### 1. Load packages and data, combine data from different tables ###
#####################################################################

# load required packages:
library(survival)
library(survminer)
library(R.matlab)
library(RColorBrewer)

# Data info:
# data from the paper by Robinson et al. (2015)
# downloaded from: https://datadryad.org/stash/dataset/doi:10.5061/dryad.m1n03

# load all data:
table1 <- read_csv("./Data/Albinama_Table1.csv")
table2_po <- read_csv("./Data/Albinama_Table2_Po.csv")
table2_pv_clin <- read_csv("./Data/Albinama_Table2_Pv_clin.csv")
table2_pv_lm <- read_csv("./Data/Albinama_Table2_Pv_LM.csv")
table2_pv_pcr <- read_csv("./Data/Albinama_Table2_Pv_PCR.csv")
table2_pvg <- read_csv("./Data/Albinama_Table2_Pvg.csv")
table3a <- read_csv("./Data/Albinama_Table3_0-3mo.csv")
table3b <- read_csv("./Data/Albinama_Table3_0-8mo.csv")
table3c <- read_csv("./Data/Albinama_Table3_4-8mo.csv")
table4_pf_clin <- read_csv("./Data/Albinama_Table4_Pf_clin.csv")
table4_pf_lm <- read_csv("./Data/Albinama_Table4_Pf_LM.csv")
table4_pf_pcr <- read_csv("./Data/Albinama_Table4_Pf_PCR.csv")
table4_pm_pcr <- read_csv("./Data/Albinama_Table4_Pm_PCR.csv")


# Combine data from different tables 

# combine data from the different tables:
data <- table1

# make age groups and add them to the data:
# age groups based on quartiles:
age_groups_qu <- ifelse(data$age<summary(data$age)[2],"4.7-6.33",ifelse(data$age<summary(data$age)[3],"6.331-7.55",ifelse(data$age<summary(data$age)[5],"7.551-8.97","8.971-10.5")))
data <- cbind(data,age_groups_qu)
# age groups equidistant:
age_groups_eq <- ifelse(data$age<6.1975,"4.7-6.1975",ifelse(data$age<7.614,"6.1975-7.614",ifelse(data$age<9.0305,"7.614-9.0305","9.0305-10.5")))
data <- cbind(data,age_groups_eq)

# P. vivax:
time_pv_pcr <- rep(NA,dim(data)[1],1)
event_pv_pcr <- rep(NA,dim(data)[1],1)
time_pv_lm <- rep(NA,dim(data)[1],1)
event_pv_lm <- rep(NA,dim(data)[1],1)
time_pv_clin <- rep(NA,dim(data)[1],1)
event_pv_clin <- rep(NA,dim(data)[1],1)
time_pv_pcr_lm <- rep(NA,dim(data)[1],1)
time_pv_pcr_clin <- rep(NA,dim(data)[1],1)
time_pv_g <- rep(NA,dim(data)[1],1)
event_pv_g <- rep(NA,dim(data)[1],1)

# P. falciparum:
time_pf_pcr <- rep(NA,dim(data)[1],1)
event_pf_pcr <- rep(NA,dim(data)[1],1)
time_pf_lm <- rep(NA,dim(data)[1],1)
event_pf_lm <- rep(NA,dim(data)[1],1)
time_pf_clin <- rep(NA,dim(data)[1],1)
event_pf_clin <- rep(NA,dim(data)[1],1)
time_pf_pcr_lm <- rep(NA,dim(data)[1],1)
time_pf_pcr_clin <- rep(NA,dim(data)[1],1)

# P. ovale:
time_po_pcr <- rep(NA,dim(data)[1],1)
event_po_pcr <- rep(NA,dim(data)[1],1)

# P. malariae:
time_pm_pcr <- rep(NA,dim(data)[1],1)
event_pm_pcr <- rep(NA,dim(data)[1],1)

for(i in 1:dim(data)[1]){
  # P. vivax: 
  if(any(table2_pv_pcr$studyid==data$studyid[i])){
    ind <- which(table2_pv_pcr$studyid==data$studyid[i])
    time_pv_pcr[i] <- table2_pv_pcr$exit_day_8mo[ind]
    event_pv_pcr[i] <- table2_pv_pcr$fail[ind]
  }else{
    time_pv_pcr[i] <- NA
    event_pv_pcr[i] <- NA
  }
  
  if(any(table2_pv_lm$studyid==data$studyid[i])){
    ind <- which(table2_pv_lm$studyid==data$studyid[i])
    time_pv_lm[i] <- table2_pv_lm$exit_day_8mo[ind]
    event_pv_lm[i] <- table2_pv_lm$fail[ind]
  }else{
    time_pv_lm[i] <- NA
    event_pv_lm[i] <- NA
  }
  
  if(any(table2_pv_clin$studyid==data$studyid[i])){
    ind <- which(table2_pv_clin$studyid==data$studyid[i])
    time_pv_clin[i] <- table2_pv_clin$exit_day_8mo[ind]
    event_pv_clin[i] <- table2_pv_clin$fail[ind]
  }else{
    time_pv_clin[i] <- NA
    event_pv_clin[i] <- NA
  }
  
  time_pv_pcr_lm[i] <- NA
  if(!is.na(event_pv_pcr[i]) & event_pv_pcr[i]==1 & !is.na(event_pv_lm[i]) & event_pv_lm[i]==1){
    time_pv_pcr_lm[i] <- time_pv_lm[i]-time_pv_pcr[i]
  }
  
  time_pv_pcr_clin[i] <- NA
  if(event_pv_pcr[i]==1 & event_pv_clin[i]==1){
    time_pv_pcr_clin[i] <- time_pv_clin[i]-time_pv_pcr[i]
  }
  
  if(any(table2_pvg$studyid==data$studyid[i])){
    ind <- which(table2_pvg$studyid==data$studyid[i])
    time_pv_g[i] <- table2_pvg$exit_day_8mo[ind]
    event_pv_g[i] <- table2_pvg$fail[ind]
  }else{
    time_pv_g[i] <- NA
    event_pv_g[i] <- NA
  }
  
  # P. falciparum: 
  if(any(table4_pf_pcr$studyid==data$studyid[i])){
    ind <- which(table4_pf_pcr$studyid==data$studyid[i])
    time_pf_pcr[i] <- table4_pf_pcr$exit_day_8mo[ind]
    event_pf_pcr[i] <- table4_pf_pcr$fail[ind]
  }else{
    time_pf_pcr[i] <- NA
    event_pf_pcr[i] <- NA
  }
  
  if(any(table4_pf_lm$studyid==data$studyid[i])){
    ind <- which(table4_pf_lm$studyid==data$studyid[i])
    time_pf_lm[i] <- table4_pf_lm$exit_day_8mo[ind]
    event_pf_lm[i] <- table4_pf_lm$fail[ind]
  }else{
    time_pf_lm[i] <- NA
    event_pf_lm[i] <- NA
  }
  
  if(any(table4_pf_clin$studyid==data$studyid[i])){
    ind <- which(table4_pf_clin$studyid==data$studyid[i])
    time_pf_clin[i] <- table4_pf_clin$exit_day_8mo[ind]
    event_pf_clin[i] <- table4_pf_clin$fail[ind]
  }else{
    time_pf_clin[i] <- NA
    event_pf_clin[i] <- NA
  }
  
  time_pf_pcr_lm[i] <- NA
  if(!is.na(event_pf_pcr[i]) & event_pf_pcr[i]==1 & !is.na(event_pf_lm[i]) & event_pf_lm[i]==1){
    time_pf_pcr_lm[i] <- time_pf_lm[i]-time_pf_pcr[i]
  }
  
  time_pf_pcr_clin[i] <- NA
  if(event_pf_pcr[i]==1 & event_pf_clin[i]==1){
    time_pf_pcr_clin[i] <- time_pf_clin[i]-time_pf_pcr[i]
  }
  
  # P. ovale:
  if(any(table2_po$studyid==data$studyid[i])){
    ind <- which(table2_po$studyid==data$studyid[i])
    time_po_pcr[i] <- table2_po$exit_day_8mo[ind]
    event_po_pcr[i] <- table2_po$fail[ind]
  }else{
    time_po_pcr[i] <- NA
    event_po_pcr[i] <- NA
  }
  
  # P. malariae:
  if(any(table4_pm_pcr$studyid==data$studyid[i])){
    ind <- which(table4_pm_pcr$studyid==data$studyid[i])
    time_pm_pcr[i] <- table4_pm_pcr$exit_day_8mo[ind]
    event_pm_pcr[i] <- table4_pm_pcr$fail[ind]
  }else{
    time_pm_pcr[i] <- NA
    event_pm_pcr[i] <- NA
  }
}

data <- cbind(data,time_pv_pcr,event_pv_pcr,time_pv_lm,event_pv_lm,time_pv_clin,event_pv_clin,time_pv_pcr_lm,time_pv_pcr_clin,time_pv_g,event_pv_g,
              time_pf_pcr,event_pf_pcr,time_pf_lm,event_pf_lm,time_pf_clin,event_pf_clin,time_pf_pcr_lm,time_pf_pcr_clin,
              time_po_pcr,event_po_pcr,time_pm_pcr,event_pm_pcr)

data$pv_pos <- as.factor(data$pv_pos)
data$pf_pos <- as.factor(data$pf_pos)
data$treat2 <- as.factor(data$treat2)
data$treat2 <- relevel(data$treat2, ref = "1")
data$villagegroup_corrected2 <- as.factor(data$villagegroup_corrected2)
data$villagegroup_corrected2 <- relevel(data$villagegroup_corrected2, ref = "2")

save(data,file='./Data/Combined_data.RData')


##########################
### 2. Survival curves ###
##########################

# all survival curves are for P. vivax infection by PCR

### Time to first recurrence for all villages combined (Fig. 1B)
table2_pv_pcr <- read_csv("./Data/Albinama_Table2_Pv_PCR.csv")
fit_pv_pcr <- survfit(Surv(exit_day_8mo,fail) ~ treat2, data = table2_pv_pcr)
ggsurvplot(fit_pv_pcr, data = table2_pv_pcr, title=expression(paste("Time to first recurrence")), pval = TRUE, 
           conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=50,xlim=c(0,230),
           legend.title = "", legend.labs = c("Placebo","Primaquine"), legend="bottom",
           ylab="Fraction without recurrence",xlab="Time [days]")

### Time to first recurrence by village (Fig. 5 confidence region)
for(i in 1:5){
  data_tmp <- table2_pv_pcr[table2_pv_pcr$villagegroup_corrected2==i,]
  fit_pv_pcr <- survfit(Surv(exit_day_8mo,fail) ~ treat2, data = data_tmp)
  summary(fit_pv_pcr)
  print(ggsurvplot(fit_pv_pcr, data = data_tmp, title=paste("Time to P. vivax infection (PCR positive, village ",i,")",sep = ""), pval = TRUE, 
                   conf.int = TRUE,risk.table = TRUE,tables.height = 0.2,tables.theme = theme_cleantable(),break.time.by=50,
                   legend.title = "", legend.labs = c("Placebo","Primaquine"), legend="bottom",
                   ylab="Fraction without P. vivax infection",xlab="Time [days]"))
}


################################
### 3. Weekly incidence rate ###
################################

### Visualize the weekly incidence rate per patients at risk (Fig. 1D)

# Use combined data:
load('./Data/Combined_data.RData')

# Incidence rate per people at risk
dt <- 7 # time steps
inc_bl <- cbind(seq(dt,max(data$time_pv_pcr[data$event_pv_pcr==1],na.rm = TRUE),dt),matrix(0,length(seq(dt,max(data$time_pv_pcr[data$event_pv_pcr==1],na.rm = TRUE),dt))))
inc_pmq <- cbind(seq(dt,max(data$time_pv_pcr[data$event_pv_pcr==1],na.rm = TRUE),dt),matrix(0,length(seq(dt,max(data$time_pv_pcr[data$event_pv_pcr==1],na.rm = TRUE),dt))))
for(i in 1:dim(inc_bl)[1]){
  if(i==1){
    inc_bl[i,2] <- sum(data$time_pv_pcr[data$treat2==0 & data$event_pv_pcr==1 & !is.na(data$time_pv_pcr)]<=inc_bl[i,1])/
      sum(data$time_pv_pcr[data$treat2==0 & !is.na(data$time_pv_pcr)]>=0) # new cases/number of patients at risk
    inc_pmq[i,2] <- sum(data$time_pv_pcr[data$treat2==1 & data$event_pv_pcr==1 & !is.na(data$time_pv_pcr)]<=inc_pmq[i,1])/
      sum(data$time_pv_pcr[data$treat2==1 & !is.na(data$time_pv_pcr)]>=0) # new cases/number of patients at risk
  }else{
    inc_bl[i,2] <- sum(data$time_pv_pcr[data$treat2==0 & data$event_pv_pcr==1 & !is.na(data$time_pv_pcr)]>=inc_bl[i-1,1] & 
                         data$time_pv_pcr[data$treat2==0 & data$event_pv_pcr==1 & !is.na(data$time_pv_pcr)]<=inc_bl[i,1])/
      sum(data$time_pv_pcr[data$treat2==0 & !is.na(data$time_pv_pcr)]>=inc_bl[i-1,1]) # new cases/number of patients at risk
    inc_pmq[i,2] <- sum(data$time_pv_pcr[data$treat2==1 & data$event_pv_pcr==1 & !is.na(data$time_pv_pcr)]>=inc_pmq[i-1,1] & 
                          data$time_pv_pcr[data$treat2==1 & data$event_pv_pcr==1 & !is.na(data$time_pv_pcr)]<=inc_pmq[i,1])/
      sum(data$time_pv_pcr[data$treat2==1 & !is.na(data$time_pv_pcr)]>=inc_pmq[i-1,1]) # new cases/number of patients at risk
  }
}

# use smooth spline instead of cubic spline:
require(splines)
data_bl <- as.data.frame(inc_bl)
data_pmq <- as.data.frame(inc_pmq)
fit_bl <- smooth.spline(data_bl$V1,data_bl$V2,df=6)
# summary(fit_bl)
fit_pmq <- smooth.spline(data_pmq$V1,data_pmq$V2,df=6)
# summary(fit_pmq)

plot(inc_bl,pch=16,col="#F8766D",las=1,xlab="Time [days]",ylab="Weekly incidence rate",
     main="Weekly incidence rate per patients at risk",xlim=c(0,max(inc_bl[,1],inc_pmq[,1])))
lines(fit_bl,col="#F8766D",lwd=2)
points(inc_pmq,col="#00BFC4",pch=16)
lines(fit_pmq,col="#00BFC4",lwd=2)
legend("topright",legend=c("primaquine", "placebo"),col=c("#00BFC4","#F8766D"), lty = 1, lwd=2)



##########################################################
### 4. Transmission and relapse rate (from model fits) ###
##########################################################

### Visualize the relationship between the relapse rate and the infection rate 
### for the 5 villages in the PNG data (Fig. 5F)

# Parameter estimates from the model fit of the population heterogeneity model 
# to all villages simultaneously with the same drug washout time distribution:
tr_v1 <- 0.002073439849889 # transmission rate
tr_v2 <- 0.000107921583196
tr_v3 <- 0.005095707640667
tr_v4 <- 0.000881235324295
tr_v5 <- 0.010107520276450

rel_v1 <- c(-5.029239751189614,1.845373879451467) # parameters for the relapse rate distribution
rel_v2 <- c(-5.694800411642519,1.632290670765690)
rel_v3 <- c(-3.574578993884039,0.935013179874915)
rel_v4 <- c(-4.393291614830225,2.258679027251658)
rel_v5 <- c(-1.235438233197847,1.127999032231824)

# Quartiles from the relapse rate distributions:
rel_q_v1 <- qlnorm(c(0.25,0.5,0.75), mean = rel_v1[1], sd = rel_v1[2])
rel_q_v2 <- qlnorm(c(0.25,0.5,0.75), mean = rel_v2[1], sd = rel_v2[2])
rel_q_v3 <- qlnorm(c(0.25,0.5,0.75), mean = rel_v3[1], sd = rel_v3[2])
rel_q_v4 <- qlnorm(c(0.25,0.5,0.75), mean = rel_v4[1], sd = rel_v4[2])
rel_q_v5 <- qlnorm(c(0.25,0.5,0.75), mean = rel_v5[1], sd = rel_v5[2])

cols <- brewer.pal(6, "Dark2") # colors

plot(c(tr_v1,tr_v2,tr_v3,tr_v4,tr_v5),c(rel_q_v1[2],rel_q_v2[2],rel_q_v3[2],rel_q_v4[2],rel_q_v5[2]),pch=16,las=1,
     xlab="Infection rate [per day]", ylab="Relapse rate [per day]",main="Relapse rate and transmission intensity for different villages",
     col = cols[2:6], ylim = c(0.001,0.63),log="y")
lines(c(tr_v1,tr_v1),rel_q_v1[c(1,3)],col=cols[2])
lines(c(tr_v2,tr_v2),rel_q_v2[c(1,3)],col=cols[3])
lines(c(tr_v3,tr_v3),rel_q_v3[c(1,3)],col=cols[4])
lines(c(tr_v4,tr_v4),rel_q_v4[c(1,3)],col=cols[5])
lines(c(tr_v5,tr_v5),rel_q_v5[c(1,3)],col=cols[6])
legend("topleft",legend=c("Village 1","Village 2","Village 3","Village 4","Village 5"),pch=16,col=cols[2:6])

# Correlation between log-transformed median relapse and transmission rate:
cor.test(log(c(rel_q_v1[2],rel_q_v2[2],rel_q_v3[2],rel_q_v4[2],rel_q_v5[2])),c(tr_v1,tr_v2,tr_v3,tr_v4,tr_v5), method="pearson")
cor.test(log(c(rel_q_v1[2],rel_q_v2[2],rel_q_v3[2],rel_q_v4[2],rel_q_v5[2])),c(tr_v1,tr_v2,tr_v3,tr_v4,tr_v5), method="spearman")

