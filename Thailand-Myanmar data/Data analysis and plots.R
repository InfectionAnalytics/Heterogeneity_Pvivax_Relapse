### Data analysis of multiple recurrence P. vivax data


### Overview:
# 1. Load packages and data
# 2. Contribution to recurrences by the number of recurrences
# 3. First and second recurrence time
# 4. Number of recurrences by time to first recurrence


#################################
### 1. Load packages and data ###
#################################

# load required packages:
library(survival)

# load data:
load("Combined_Time_Event.RData")

# save data in a MATLAB file
# library(R.matlab)
# writeMat("Combined_Time_Data.mat",Combined_Time_Data=Combined_Time_Data)


###################################################################
### 2. Contribution to recurrences by the number of recurrences ###
###################################################################


# consider all data:
data <- Combined_Time_Data # For Fig. S6
# consider only the data for individuals with blood-stage treatment only (no primaquine):
# data <- Combined_Time_Data[which(Combined_Time_Data$arm_num!="CHQ/PMQ"),]

# number of recurrences for each individual:
tmp <- unique(data[data$Censored==0,c(1,2)])
ind_num_rec <- sapply(c(1:length(unique(tmp$patientid))),function(i){max(tmp$episode[which(tmp$patientid==unique(tmp$patientid)[i])])-1})#[unique(tmp$patientid)[i],max(tmp$episode[which(tmp$patientid==unique(tmp$patientid)[i])])-1]})
ind_num_rec <- as.data.frame(cbind(unique(tmp$patientid),ind_num_rec))
colnames(ind_num_rec) <- c("id","num_rec")
ind_num_rec$num_rec <- as.numeric(as.character(ind_num_rec$num_rec))

# number and percent of individuals by number of recurrences:
n_ind_by_rec <- sapply(c(1:14),function(i){length(ind_num_rec$id[which(ind_num_rec$num_rec==i)])})
n_ind_by_rec <- c(length(unique(data$patientid))-length(ind_num_rec$id),n_ind_by_rec)
perc_ind_by_rec <- n_ind_by_rec*100/sum(n_ind_by_rec)

# number and percent of recurrences caused by individuals with a certain number of recurrences:
n_rec_by_rec <- rep(0,max(ind_num_rec$num_rec))
for(i in 1:max(ind_num_rec$num_rec)){
  tmp_id <- ind_num_rec$id[which(ind_num_rec$num_rec==i)]
  n_rec_by_rec[i] <- dim(data[which(data$patientid%in%tmp_id & data$episode>1 & data$Censored==0),])[1]
}
n_rec_by_rec <- c(0,n_rec_by_rec)
perc_rec_by_rec <- n_rec_by_rec*100/sum(n_rec_by_rec)

# combined results (Table S10):
table_contributions <- data.frame(n_rec=c(0:14),n_ind=n_ind_by_rec,perc_ind=perc_ind_by_rec,n_rec_by_rec,perc_rec_by_rec)

# cumulative distributions (from highest to lowest number of recurrences):
cum_ind_by_rec <- sapply(c(1:length(perc_ind_by_rec)),function(i){sum(perc_ind_by_rec[(length(perc_ind_by_rec)-i):length(perc_ind_by_rec)])})/100
cum_rec_by_rec <- sapply(c(1:length(perc_rec_by_rec)),function(i){sum(perc_rec_by_rec[(length(perc_rec_by_rec)-i):length(perc_rec_by_rec)])})/100

# plot contribution to recurrences
plot(cum_ind_by_rec,cum_rec_by_rec,type="l",lwd=2,las=1,xlab="Proportion of the population",ylab="Proportion of recurrences",main="Contribution to recurrences")
points(cum_ind_by_rec,cum_rec_by_rec,pch=16)


###########################################
### 3. First and second recurrence time ###
###########################################


# Time to first and second recurrence for all individuals with at least one recurrence
data_tmp <- cbind(Combined_Time_Data,study=substr(Combined_Time_Data$patientid,1,3))
data_tmp$arm_num[data_tmp$PMQ_partner=="DP"] <- "DP/PMQ"
data_tmp <- cbind(data_tmp,drug_study=paste(data_tmp$arm_num,"_",data_tmp$study,sep=""))
id <- data_tmp$patientid[which(data_tmp$episode==3)] # individuals with at least one recurrence
time1 <- data_tmp$Time_to_event[data_tmp$patientid%in%id & data_tmp$episode==2] # time to first recurrence
event1 <- 1-data_tmp$Censored[data_tmp$patientid%in%id & data_tmp$episode==2] # event at first recurrence
time2 <- data_tmp$Time_to_event[data_tmp$patientid%in%id & data_tmp$episode==3] # time to second recurrence
event2 <- 1-data_tmp$Censored[data_tmp$patientid%in%id & data_tmp$episode==3] # event at second recurrence
drug <- data_tmp$arm_num[data_tmp$patientid%in%id & data_tmp$episode==2]
study <- data_tmp$study[data_tmp$patientid%in%id & data_tmp$episode==2]
drug_study <- data_tmp$drug_study[data_tmp$patientid%in%id & data_tmp$episode==2]
data_tmp <- data.frame(id=id,drug=drug,time1=as.numeric(time1),event1=as.numeric(event1),time2=as.numeric(time2),event2=as.numeric(event2),study,drug_study)
data_tmp_vhx <- data_tmp[data_tmp$study=="VHX",]
data_tmp_bpd <- data_tmp[data_tmp$study=="BPD",]
data_tmp_as <- data_tmp[data_tmp$drug=="AS",]
data_tmp_chq <- data_tmp[data_tmp$drug=="CHQ",]
data_tmp_pmq <- data_tmp[data_tmp$drug%in%c("CHQ/PMQ","DP/PMQ"),]


# Correlations between first and second recurrence time
# with censored second recurrences
corel_all <- cor(data_tmp$time1,data_tmp$time2,method="spearman")
corel_all_pval <- cor.test(data_tmp$time1,data_tmp$time2,method = "spearman")$p.value
corel_vhx <- cor(data_tmp_vhx$time1,data_tmp_vhx$time2,method="spearman")
corel_vhx_pval <- cor.test(data_tmp_vhx$time1,data_tmp_vhx$time2,method = "spearman")$p.value
corel_bpd <- cor(data_tmp_bpd$time1,data_tmp_bpd$time2,method="spearman")
corel_bpd_pval <- cor.test(data_tmp_bpd$time1,data_tmp_bpd$time2,method = "spearman")$p.value
corel_as <- cor(data_tmp_as$time1,data_tmp_as$time2,method="spearman")
corel_as_pval <- cor.test(data_tmp_as$time1,data_tmp_as$time2,method = "spearman")$p.value
corel_chq <- cor(data_tmp_chq$time1,data_tmp_chq$time2,method="spearman")
corel_chq_pval <- cor.test(data_tmp_chq$time1,data_tmp_chq$time2,method = "spearman")$p.value
corel_pmq <- cor(data_tmp_pmq$time1,data_tmp_pmq$time2,method="spearman")
corel_pmq_pval <- cor.test(data_tmp_pmq$time1,data_tmp_pmq$time2,method = "spearman")$p.value

table_cor <- data.frame(data=c("all","VHX","BPD","AS","CHQ","PMQ+"),corr_all=c(corel_all,corel_vhx,corel_bpd,corel_as,corel_chq,corel_pmq),
                        corr_all_p=c(corel_all_pval,corel_vhx_pval,corel_bpd_pval,corel_as_pval,corel_chq_pval,corel_pmq_pval))

# without censored second recurrences
corel_all <- cor(data_tmp$time1[data_tmp$event2==1],data_tmp$time2[data_tmp$event2==1],method="spearman")
corel_all_pval <- cor.test(data_tmp$time1[data_tmp$event2==1],data_tmp$time2[data_tmp$event2==1],method = "spearman")$p.value
corel_vhx <- cor(data_tmp_vhx$time1[data_tmp_vhx$event2==1],data_tmp_vhx$time2[data_tmp_vhx$event2==1],method="spearman")
corel_vhx_pval <- cor.test(data_tmp_vhx$time1[data_tmp_vhx$event2==1],data_tmp_vhx$time2[data_tmp_vhx$event2==1],method = "spearman")$p.value
corel_bpd <- cor(data_tmp_bpd$time1[data_tmp_bpd$event2==1],data_tmp_bpd$time2[data_tmp_bpd$event2==1],method="spearman")
corel_bpd_pval <- cor.test(data_tmp_bpd$time1[data_tmp_bpd$event2==1],data_tmp_bpd$time2[data_tmp_bpd$event2==1],method = "spearman")$p.value
corel_as <- cor(data_tmp_as$time1[data_tmp_as$event2==1],data_tmp_as$time2[data_tmp_as$event2==1],method="spearman")
corel_as_pval <- cor.test(data_tmp_as$time1[data_tmp_as$event2==1],data_tmp_as$time2[data_tmp_as$event2==1],method = "spearman")$p.value
corel_chq <- cor(data_tmp_chq$time1[data_tmp_chq$event2==1],data_tmp_chq$time2[data_tmp_chq$event2==1],method="spearman")
corel_chq_pval <- cor.test(data_tmp_chq$time1[data_tmp_chq$event2==1],data_tmp_chq$time2[data_tmp_chq$event2==1],method = "spearman")$p.value
corel_pmq <- cor(data_tmp_pmq$time1[data_tmp_pmq$event2==1],data_tmp_pmq$time2[data_tmp_pmq$event2==1],method="spearman")
corel_pmq_pval <- cor.test(data_tmp_pmq$time1[data_tmp_pmq$event2==1],data_tmp_pmq$time2[data_tmp_pmq$event2==1],method = "spearman")$p.value

# combine everything (Table S11 & Table 1 for the Thailand-Myanmar data):
table_cor <- cbind(table_cor,corr_no_cens=c(corel_all,corel_vhx,corel_bpd,corel_as,corel_chq,corel_pmq),
                   corr_no_cens_p=c(corel_all_pval,corel_vhx_pval,corel_bpd_pval,corel_as_pval,corel_chq_pval,corel_pmq_pval))

# Plot second recurrence time vs first recurrence time by study: (Fig. S5A)
plot(data_tmp_vhx$time1[which(data_tmp_vhx$event2==1)],data_tmp_vhx$time2[which(data_tmp_vhx$event2==1)],pch=16,col="black",cex=1,
     main="Association between time to first recurrence and time between first and second recurrence",xlab="Time to first recurrence [days]",
     ylab="Time between first and second recurrence [days]",xlim=c(0,370),ylim = c(0,370))
points(data_tmp_vhx$time1[which(data_tmp_vhx$event2==0)],data_tmp_vhx$time2[which(data_tmp_vhx$event2==0)],pch=16,col="red",cex=1)
points(data_tmp_bpd$time1[which(data_tmp_bpd$event2==1)],data_tmp_bpd$time2[which(data_tmp_bpd$event2==1)],pch=17,col="black",cex=1)
points(data_tmp_bpd$time1[which(data_tmp_bpd$event2==0)],data_tmp_bpd$time2[which(data_tmp_bpd$event2==0)],pch=17,col="red",cex=1)
legend("topright", legend=c("VHX study", "VHX study, right censored", "BPD study", "BPD study, right censored"), 
       col=c("black", "red", "black", "red"), pch=c(16,16,17,17), cex=1)

# Plot second recurrence time vs first recurrence time by drug: (Fig. S5B)
plot(data_tmp_as$time1[which(data_tmp_as$event2==1)],data_tmp_as$time2[which(data_tmp_as$event2==1)],pch=16,col="black",cex=1,
     main="Association between time to first recurrence and time between first and second recurrence",xlab="Time to first recurrence [days]",
     ylab="Time between first and second recurrence [days]",xlim=c(0,370),ylim = c(0,370),las=1)
points(data_tmp_as$time1[which(data_tmp_as$event2==0)],data_tmp_as$time2[which(data_tmp_as$event2==0)],pch=16,col="red",cex=1)
points(data_tmp_chq$time1[which(data_tmp_chq$event2==1)],data_tmp_chq$time2[which(data_tmp_chq$event2==1)],pch=17,col="black",cex=1)
points(data_tmp_chq$time1[which(data_tmp_chq$event2==0)],data_tmp_chq$time2[which(data_tmp_chq$event2==0)],pch=17,col="red",cex=1)
points(data_tmp_pmq$time1[which(data_tmp_pmq$event2==1)],data_tmp_pmq$time2[which(data_tmp_pmq$event2==1)],pch=18,col="black",cex=1)
points(data_tmp_pmq$time1[which(data_tmp_pmq$event2==0)],data_tmp_pmq$time2[which(data_tmp_pmq$event2==0)],pch=18,col="red",cex=1)
legend("topright", legend=c("Artesunate", "Artesunate, right censored", "Chloroquine", "Chloroquine, right censored","Primaquine+",
                            "Primaquine+, right censored"), col=c("black", "red", "black", "red","black","red"), pch=c(16,16,17,17,18,18), cex=1)


# Correlations with recurrence times limited (Thailand-Myanmar data in Table S12):
time_limit <- 182 # time limit 182 days
ind_as <- which(data_tmp_as$time1 <= time_limit & data_tmp_as$time2<=time_limit & data_tmp_as$event2==1)
corel_as_half_year <- cor(data_tmp_as$time1[ind_as],data_tmp_as$time2[ind_as],method="spearman")
corel_as_half_year_pval <- cor.test(data_tmp_as$time1[ind_as],data_tmp_as$time2[ind_as],method = "spearman")$p.value
ind_chq <- which(data_tmp_chq$time1 <= time_limit & data_tmp_chq$time2<=time_limit & data_tmp_chq$event2==1)
corel_chq_half_year <- cor(data_tmp_chq$time1[ind_chq],data_tmp_chq$time2[ind_chq],method="spearman")
corel_chq_half_year_pval <- cor.test(data_tmp_chq$time1[ind_chq],data_tmp_chq$time2[ind_chq],method = "spearman")$p.value


# Cox regression on time from first to second recurrence (Table S13):
fit.coxph_as <- coxph(Surv(time2,event2) ~ time1, data = data_tmp_as)
summary(fit.coxph_as)
exp(confint(fit.coxph_as))

fit.coxph_chq <- coxph(Surv(time2,event2) ~ time1, data = data_tmp_chq)
summary(fit.coxph_chq)
exp(confint(fit.coxph_chq))


# Likelihood-ratio test for comparing model fits 3 & 4 to the first and second recurrence time:
nllh.3 <- 4203 # negative log-likelihood for model 3
n.par.3 <- 10 # number of parameters for model 3
nllh.4 <- 4191 # negative log-likelihood for model 4
n.par.4 <- 11 # number of parameters for model 4

lr <- 2*((-nllh.4)-(-nllh.3))
p.val <- pchisq(lr, df = n.par.4-n.par.3, lower.tail = FALSE)


############################################################
### 4. Number of recurrences by time to first recurrence ###
############################################################


# Number of recurrences per year by time to first recurrence for each individual with follow-up of at least one year
# use data grouped by drug and study, only use AS and CHQ treated individuals for this data analysis
data_tmp <- Combined_Time_Data[Combined_Time_Data$arm_num%in%c("AS","CHQ"),]
# patientid, time to first recurrence and number of recurrences within 1 year for each individual
n_rec_AS <- matrix(NA,sum(data_tmp$Censored[data_tmp$arm_num=="AS" & data_tmp$episode==2 & data_tmp$FU_time>=365]==0),3)
n_rec_CHQ <- matrix(NA,sum(data_tmp$Censored[data_tmp$arm_num=="CHQ" & data_tmp$episode==2 & data_tmp$FU_time>=365]==0),3)

for(i in 1:length(unique(data_tmp$patientid))){ # for each patient
  ind <- unique(data_tmp$patientid)[i]
  if(any(data_tmp$patientid==ind & data_tmp$episode==2 & data_tmp$Censored==0 & data_tmp$FU_time>=365)){# if there is a first event
    drug <- unique(data_tmp$arm_num[data_tmp$patientid==ind])
    time_1st <- data_tmp$Time_to_event[data_tmp$patientid==ind & data_tmp$episode==2 & data_tmp$Censored==0]
    n_rec <- sum(data_tmp$patientid==ind & data_tmp$Censored==0 & data_tmp$Time_since_enrolment<=365) # number of all recurrences during the first 365 days of follow-up
    if(drug=="AS"){
      n_rec_AS[which(is.na(n_rec_AS[,1]))[1],] <- c(ind,time_1st,n_rec)
    }else if(drug=="CHQ"){
      n_rec_CHQ[which(is.na(n_rec_CHQ[,1]))[1],] <- c(ind,time_1st,n_rec)
    }
  }
}

plot(n_rec_AS[,2],n_rec_AS[,3],pch=16,main="AS",xlab="Time to first recurrence [days]",ylab="Number of recurrences",las=1)
plot(n_rec_CHQ[,2],n_rec_CHQ[,3],pch=16,main="CHQ",xlab="Time to first recurrence [days]",ylab="Number of recurrences",las=1)

# fit a smooth spline through the data:
require(splines)
data_as <- as.data.frame(n_rec_AS)
data_as$V2 <- as.numeric(as.character(data_as$V2))
data_as$V3 <- as.numeric(as.character(data_as$V3))
data_chq <- as.data.frame(n_rec_CHQ)
data_chq$V2 <- as.numeric(as.character(data_chq$V2))
data_chq$V3 <- as.numeric(as.character(data_chq$V3))
fit_as<-smooth.spline(data_as$V2,data_as$V3,df=4)
fit_chq <-smooth.spline(data_chq$V2,data_chq$V3,df=4)

# Fig. S13A for the data:
plot(n_rec_AS[,2],n_rec_AS[,3],pch=16,main="AS",xlab="Time to first recurrence [days]",ylab="Number of recurrences",las=1,xlim=c(0,365),ylim=c(1,17))
lines(fit_as,col="black",lwd=2)

# Fig. S13B for the data:
plot(n_rec_CHQ[,2],n_rec_CHQ[,3],pch=16,main="CHQ",xlab="Time to first recurrence [days]",ylab="Number of recurrences",las=1,xlim=c(0,365))
lines(fit_chq,col="black",lwd=2)





