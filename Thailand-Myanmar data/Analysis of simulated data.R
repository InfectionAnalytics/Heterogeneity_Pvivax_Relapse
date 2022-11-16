### Analysis of simulated multiple recurrence data


### Overview:
# 1. Load packages and data
# 2. Correlation between first and second recurrence times
# 3. Survival curves
# 4. Average number of recurrences per individual and year by time to first recurrence


#################################
### 1. Load packages and data ###
#################################

# load required packages:
library(survival)
library(survminer)
library(R.matlab)

### load simulated data from MATLAB
# Models:
# 1: constant relapse rate
# 2: temporal heterogeneity
# 3: population heterogeneity
data_AS <- readMat('./Results/Simulation_1year_AS_VHX (1000 sim, mod 1 to 3).mat')
sim1_AS <- data_AS$sim1
sim2_AS <- data_AS$sim2
sim3_AS <- data_AS$sim3

data_CHQ <- readMat('./Results/Simulation_1year_CHQ_VHX (1000 sim, mod 1 to 3).mat')
sim1_CHQ <- data_CHQ$sim1
sim2_CHQ <- data_CHQ$sim2
sim3_CHQ <- data_CHQ$sim3


################################################################
### 2. Correlation between first and second recurrence times ###
################################################################

# save first and second recurrence times (events only)
times_AS_1  <- NULL
times_AS_2  <- NULL
times_AS_3  <- NULL
times_CHQ_1  <- NULL
times_CHQ_2  <- NULL
times_CHQ_3  <- NULL

for(l in 1:1000){
  sim1_AS_tmp <- sim1_AS[[l]][[1]]
  ind <- unique(sim1_AS_tmp[which(sim1_AS_tmp[,4]==2 & sim1_AS_tmp[,2]>0),1])
  time1 <- sim1_AS_tmp[which(sim1_AS_tmp[,1] %in% ind & sim1_AS_tmp[,4]==1),2]
  time2 <- sim1_AS_tmp[which(sim1_AS_tmp[,1] %in% ind & sim1_AS_tmp[,4]==2),2]
  times_AS_1 <- rbind(times_AS_1,cbind(time1,time2))
  
  sim2_AS_tmp <- sim2_AS[[l]][[1]]
  ind <- unique(sim2_AS_tmp[which(sim2_AS_tmp[,4]==2 & sim2_AS_tmp[,2]>0),1])
  time1 <- sim2_AS_tmp[which(sim2_AS_tmp[,1] %in% ind & sim2_AS_tmp[,4]==1),2]
  time2 <- sim2_AS_tmp[which(sim2_AS_tmp[,1] %in% ind & sim2_AS_tmp[,4]==2),2]
  times_AS_2 <- rbind(times_AS_2,cbind(time1,time2))
  
  sim3_AS_tmp <- sim3_AS[[l]][[1]]
  ind <- unique(sim3_AS_tmp[which(sim3_AS_tmp[,4]==2 & sim3_AS_tmp[,2]>0),1])
  time1 <- sim3_AS_tmp[which(sim3_AS_tmp[,1] %in% ind & sim3_AS_tmp[,4]==1),2]
  time2 <- sim3_AS_tmp[which(sim3_AS_tmp[,1] %in% ind & sim3_AS_tmp[,4]==2),2]
  times_AS_3 <- rbind(times_AS_3,cbind(time1,time2))
  
  sim1_CHQ_tmp <- sim1_CHQ[[l]][[1]]
  ind <- unique(sim1_CHQ_tmp[which(sim1_CHQ_tmp[,4]==2 & sim1_CHQ_tmp[,2]>0),1])
  time1 <- sim1_CHQ_tmp[which(sim1_CHQ_tmp[,1] %in% ind & sim1_CHQ_tmp[,4]==1),2]
  time2 <- sim1_CHQ_tmp[which(sim1_CHQ_tmp[,1] %in% ind & sim1_CHQ_tmp[,4]==2),2]
  times_CHQ_1 <- rbind(times_CHQ_1,cbind(time1,time2))
  
  sim2_CHQ_tmp <- sim2_CHQ[[l]][[1]]
  ind <- unique(sim2_CHQ_tmp[which(sim2_CHQ_tmp[,4]==2 & sim2_CHQ_tmp[,2]>0),1])
  time1 <- sim2_CHQ_tmp[which(sim2_CHQ_tmp[,1] %in% ind & sim2_CHQ_tmp[,4]==1),2]
  time2 <- sim2_CHQ_tmp[which(sim2_CHQ_tmp[,1] %in% ind & sim2_CHQ_tmp[,4]==2),2]
  times_CHQ_2 <- rbind(times_CHQ_2,cbind(time1,time2))
  
  sim3_CHQ_tmp <- sim3_CHQ[[l]][[1]]
  ind <- unique(sim3_CHQ_tmp[which(sim3_CHQ_tmp[,4]==2 & sim3_CHQ_tmp[,2]>0),1])
  time1 <- sim3_CHQ_tmp[which(sim3_CHQ_tmp[,1] %in% ind & sim3_CHQ_tmp[,4]==1),2]
  time2 <- sim3_CHQ_tmp[which(sim3_CHQ_tmp[,1] %in% ind & sim3_CHQ_tmp[,4]==2),2]
  times_CHQ_3 <- rbind(times_CHQ_3,cbind(time1,time2))
}

tmp <- cor.test(times_AS_1[,1],times_AS_1[,2],method = "spearman")
c_AS_1 <- c(tmp[[4]],tmp[[3]])
tmp <- cor.test(times_AS_2[,1],times_AS_2[,2],method = "spearman")
c_AS_2 <- c(tmp[[4]],tmp[[3]])
tmp <- cor.test(times_AS_3[,1],times_AS_3[,2],method = "spearman")
c_AS_3 <- c(tmp[[4]],tmp[[3]])

tmp <- cor.test(times_CHQ_1[,1],times_CHQ_1[,2],method = "spearman")
c_CHQ_1 <- c(tmp[[4]],tmp[[3]])
tmp <- cor.test(times_CHQ_2[,1],times_CHQ_2[,2],method = "spearman")
c_CHQ_2 <- c(tmp[[4]],tmp[[3]])
tmp <- cor.test(times_CHQ_3[,1],times_CHQ_3[,2],method = "spearman")
c_CHQ_3 <- c(tmp[[4]],tmp[[3]])

# Table 1 (for models):
table.cor <- data.frame(model=c(1,2,3),as.cor=c(c_AS_1[1],c_AS_2[1],c_AS_3[1]),as.cor.p=c(c_AS_1[2],c_AS_2[2],c_AS_3[2]),
                        chq.cor=c(c_CHQ_1[1],c_CHQ_2[1],c_CHQ_3[1]),chq.cor.p=c(c_CHQ_1[2],c_CHQ_2[2],c_CHQ_3[2]))

# Correlation for recurrence times <= 182 days
time_limit <- 182 # time limit in days
ind <- which(times_AS_1[,1]<=time_limit & times_AS_1[,2]<=time_limit)
tmp <- cor.test(times_AS_1[ind,1],times_AS_1[ind,2],method="spearman")
c_AS_1_half_year <- c(tmp[[4]],tmp[[3]])

ind <- which(times_AS_2[,1]<=time_limit & times_AS_2[,2]<=time_limit)
tmp <- cor.test(times_AS_2[ind,1],times_AS_2[ind,2],method="spearman")
c_AS_2_half_year <- c(tmp[[4]],tmp[[3]])

ind <- which(times_AS_3[,1]<=time_limit & times_AS_3[,2]<=time_limit)
tmp <- cor.test(times_AS_3[ind,1],times_AS_3[ind,2],method="spearman")
c_AS_3_half_year <- c(tmp[[4]],tmp[[3]])

ind <- which(times_CHQ_1[,1]<=time_limit & times_CHQ_1[,2]<=time_limit)
tmp <- cor.test(times_CHQ_1[ind,1],times_CHQ_1[ind,2],method="spearman")
c_CHQ_1_half_year <- c(tmp[[4]],tmp[[3]])

ind <- which(times_CHQ_2[,1]<=time_limit & times_CHQ_2[,2]<=time_limit)
tmp <- cor.test(times_CHQ_2[ind,1],times_CHQ_2[ind,2],method="spearman")
c_CHQ_2_half_year <- c(tmp[[4]],tmp[[3]])

ind <- which(times_CHQ_3[,1]<=time_limit & times_CHQ_3[,2]<=time_limit)
tmp <- cor.test(times_CHQ_3[ind,1],times_CHQ_3[ind,2],method="spearman")
c_CHQ_3_half_year <- c(tmp[[4]],tmp[[3]])

# Table S12 (for models):
table.cor.half.year <- data.frame(model=c(1,2,3),as.cor=c(c_AS_1_half_year[1],c_AS_2_half_year[1],c_AS_3_half_year[1]),
                        as.cor.p=c(c_AS_1_half_year[2],c_AS_2_half_year[2],c_AS_3_half_year[2]),
                        chq.cor=c(c_CHQ_1_half_year[1],c_CHQ_2_half_year[1],c_CHQ_3_half_year[1]),
                        chq.cor.p=c(c_CHQ_1_half_year[2],c_CHQ_2_half_year[2],c_CHQ_3_half_year[2]))


##########################
### 3. Survival curves ###
##########################

# make the survival curves and save them to plot all survival curves together
# surv_curve_data_AS_mod1 <- vector("list", 1000) 
# surv_curve_data_AS_mod2 <- vector("list", 1000) 
# surv_curve_data_AS_mod3 <- vector("list", 1000) 
# surv_curve_data_CHQ_mod1 <- vector("list", 1000) 
# surv_curve_data_CHQ_mod2 <- vector("list", 1000) 
# surv_curve_data_CHQ_mod3 <- vector("list", 1000) 
# surv_curve_data_AS_c_mod1 <- vector("list", 1000) # censoring data
# surv_curve_data_AS_c_mod2 <- vector("list", 1000) 
# surv_curve_data_AS_c_mod3 <- vector("list", 1000) 
# surv_curve_data_CHQ_c_mod1 <- vector("list", 1000) 
# surv_curve_data_CHQ_c_mod2 <- vector("list", 1000) 
# surv_curve_data_CHQ_c_mod3 <- vector("list", 1000) 
# for(i in 1:1000){
#   for(m in 1:3){
#     if(m==1){
#       data_tmp_AS <- sim1_AS[[i]][[1]]
#       data_tmp_CHQ <- sim1_CHQ[[i]][[1]]
#     }else if(m==2){
#       data_tmp_AS <- sim2_AS[[i]][[1]]
#       data_tmp_CHQ <- sim2_CHQ[[i]][[1]]
#     }else if(m==3){
#       data_tmp_AS <- sim3_AS[[i]][[1]]
#       data_tmp_CHQ <- sim3_CHQ[[i]][[1]]
#     }
#     
#     event1_AS <- matrix(1,1000,1)
#     event1_AS[data_tmp_AS[which(data_tmp_AS[,4]==0),1]] <- 0
#     time1_AS <- data_tmp_AS[which(data_tmp_AS[,4]==1 | data_tmp_AS[,4]==0),2]
#     time1_AS[data_tmp_AS[which(data_tmp_AS[,4]==0),1]] <- 365
#     
#     event1_CHQ <- matrix(1,1000,1)
#     event1_CHQ[data_tmp_CHQ[which(data_tmp_CHQ[,4]==0),1]] <- 0
#     time1_CHQ <- data_tmp_CHQ[which(data_tmp_CHQ[,4]==1 | data_tmp_CHQ[,4]==0),2]
#     time1_CHQ[data_tmp_CHQ[which(data_tmp_CHQ[,4]==0),1]] <- 365
#     
#     event2_AS <- matrix(1,1000,1)
#     time2_AS <- matrix(1,1000,1)
#     event2_CHQ <- matrix(1,1000,1)
#     time2_CHQ <- matrix(1,1000,1)
#     for(j in 1:1000){
#       if(max(data_tmp_AS[data_tmp_AS[,1]==j,4])>1){
#         time2_AS[j] <- data_tmp_AS[which(data_tmp_AS[,1]==j & data_tmp_AS[,4]==2),2]
#       }else if(max(data_tmp_AS[data_tmp_AS[,1]==j,4])==1){
#         event2_AS[j] <- 0
#         time2_AS[j] <- 365-data_tmp_AS[which(data_tmp_AS[,1]==j & data_tmp_AS[,4]==1),2]
#       }else{
#         event2_AS[j] <- 0
#         time2_AS[j] <- NaN
#       }
#       if(max(data_tmp_CHQ[data_tmp_CHQ[,1]==j,4])>1){
#         time2_CHQ[j] <- data_tmp_CHQ[which(data_tmp_CHQ[,1]==j & data_tmp_CHQ[,4]==2),2]
#       }else if(max(data_tmp_CHQ[data_tmp_CHQ[,1]==j,4])==1){
#         event2_CHQ[j] <- 0
#         time2_CHQ[j] <- 365-data_tmp_CHQ[which(data_tmp_CHQ[,1]==j & data_tmp_CHQ[,4]==1),2]
#       }else{
#         event2_CHQ[j] <- 0
#         time2_CHQ[j] <- NaN
#       }
#     }
#     
#     # remove individuals with no recurrences
#     if(any(event1_AS==0)){
#       time1_AS <- time1_AS[-which(event1_AS==0)]
#       time2_AS <- time2_AS[-which(event1_AS==0)]
#       event2_AS <- event2_AS[-which(event1_AS==0)]
#     }
#     if(any(event1_CHQ==0)){
#       time1_CHQ <- time1_CHQ[-which(event1_CHQ==0)]
#       time2_CHQ <- time2_CHQ[-which(event1_CHQ==0)]
#       event2_CHQ <- event2_CHQ[-which(event1_CHQ==0)]
#     }
#     
#     # grouping of first recurrence times by quartiles:
#     groups_AS <- sapply(c(1:length(time1_AS)),function(i){ifelse(time1_AS[i]<summary(time1_AS)[2],1,ifelse(time1_AS[i]<summary(time1_AS)[3],2,
#                                                                                                            ifelse(time1_AS[i]<summary(time1_AS)[5],3,4)))})
#     groups_AS_names <- c("quartile 1","quartile 2","quartile 3","quartile 4")
#     groups_AS <- groups_AS_names[groups_AS]
#     
#     groups_CHQ <- sapply(c(1:length(time1_CHQ)),function(i){ifelse(time1_CHQ[i]<summary(time1_CHQ)[2],1,ifelse(time1_CHQ[i]<summary(time1_CHQ)[3],2,
#                                                                                                                ifelse(time1_CHQ[i]<summary(time1_CHQ)[5],3,4)))})
#     groups_CHQ_names <- c("quartile 1","quartile 2","quartile 3","quartile 4")
#     groups_CHQ <- groups_CHQ_names[groups_CHQ]
#     
#     # data for survival curves:
#     AS_surv <- as.data.frame(cbind(time2_AS,groups_AS,event2_AS))
#     names(AS_surv) <- c("time","group","event")
#     AS_surv$time <- as.numeric(AS_surv$time)
#     AS_surv$event <- as.numeric(AS_surv$event)
#     
#     CHQ_surv <- as.data.frame(cbind(time2_CHQ,groups_CHQ,event2_CHQ))
#     names(CHQ_surv) <- c("time","group","event")
#     CHQ_surv$time <- as.numeric(CHQ_surv$time)
#     CHQ_surv$event <- as.numeric(CHQ_surv$event)
#     
#     # survival curves:
#     fit_AS <- survfit(Surv(time=time,event=event) ~ group, data = AS_surv)
#     fit_CHQ <- survfit(Surv(time=time,event=event) ~ group, data = CHQ_surv)
#     
#     # save survival curve data and manually plot all curves together later
#     surv_curve_data_AS_tmp <- as.data.frame(cbind(summary(fit_AS)$time,summary(fit_AS)$surv,substr(summary(fit_AS)$strata,7,16)))
#     names(surv_curve_data_AS_tmp) <- c("time","surv","quartile")
#     surv_curve_data_AS_tmp$quartile <- factor(surv_curve_data_AS_tmp$quartile)
#     surv_curve_data_AS_tmp$time <- as.numeric(as.character(surv_curve_data_AS_tmp$time))
#     surv_curve_data_AS_tmp$surv <- as.numeric(as.character(surv_curve_data_AS_tmp$surv))
#     
#     surv_curve_data_CHQ_tmp <- as.data.frame(cbind(summary(fit_CHQ)$time,summary(fit_CHQ)$surv,substr(summary(fit_CHQ)$strata,7,16)))
#     names(surv_curve_data_CHQ_tmp) <- c("time","surv","quartile")
#     surv_curve_data_CHQ_tmp$quartile <- factor(surv_curve_data_CHQ_tmp$quartile)
#     surv_curve_data_CHQ_tmp$time <- as.numeric(as.character(surv_curve_data_CHQ_tmp$time))
#     surv_curve_data_CHQ_tmp$surv <- as.numeric(as.character(surv_curve_data_CHQ_tmp$surv))
#     
#     # save censoring data for plots
#     surv_curve_data_AS_c_tmp <- dplyr::filter(AS_surv,event==0)
#     surv <- sapply(c(1:dim(surv_curve_data_AS_c_tmp)[1]),function(x){ifelse(length(which(surv_curve_data_AS_tmp$quartile==surv_curve_data_AS_c_tmp$group[x] & surv_curve_data_AS_tmp$time<=as.numeric(surv_curve_data_AS_c_tmp$time[x])))>=1,
#                                                                             surv_curve_data_AS_tmp[max(which(surv_curve_data_AS_tmp$quartile==surv_curve_data_AS_c_tmp$group[x] & 
#                                                                                                                surv_curve_data_AS_tmp$time<=as.numeric(surv_curve_data_AS_c_tmp$time[x]))),2],1)})
#     surv_curve_data_AS_c_tmp <- as.data.frame(cbind(dplyr::select(surv_curve_data_AS_c_tmp,time,group),surv))
#     names(surv_curve_data_AS_c_tmp) <- c("time","quartile","surv")
#     surv_curve_data_AS_c_tmp$time <- as.numeric(as.character(surv_curve_data_AS_c_tmp$time))
#     surv_curve_data_AS_c_tmp$surv <- as.numeric(as.character(surv_curve_data_AS_c_tmp$surv))
#     
#     surv_curve_data_CHQ_c_tmp <- dplyr::filter(CHQ_surv,event==0)
#     surv <- sapply(c(1:dim(surv_curve_data_CHQ_c_tmp)[1]),function(x){ifelse(length(which(surv_curve_data_CHQ_tmp$quartile==surv_curve_data_CHQ_c_tmp$group[x] & surv_curve_data_CHQ_tmp$time<=as.numeric(surv_curve_data_CHQ_c_tmp$time[x])))>=1,
#                                                                              surv_curve_data_CHQ_tmp[max(which(surv_curve_data_CHQ_tmp$quartile==surv_curve_data_CHQ_c_tmp$group[x] & 
#                                                                                                                  surv_curve_data_CHQ_tmp$time<=as.numeric(surv_curve_data_CHQ_c_tmp$time[x]))),2],1)})
#     surv_curve_data_CHQ_c_tmp <- as.data.frame(cbind(dplyr::select(surv_curve_data_CHQ_c_tmp,time,group),surv))
#     names(surv_curve_data_CHQ_c_tmp) <- c("time","quartile","surv")
#     surv_curve_data_CHQ_c_tmp$time <- as.numeric(as.character(surv_curve_data_CHQ_c_tmp$time))
#     surv_curve_data_CHQ_c_tmp$surv <- as.numeric(as.character(surv_curve_data_CHQ_c_tmp$surv))
#     
#     if(m==1){
#       surv_curve_data_AS_mod1[[i]] <- surv_curve_data_AS_tmp
#       surv_curve_data_CHQ_mod1[[i]] <- surv_curve_data_CHQ_tmp
#       surv_curve_data_AS_c_mod1[[i]] <- surv_curve_data_AS_c_tmp
#       surv_curve_data_CHQ_c_mod1[[i]] <- surv_curve_data_CHQ_c_tmp
#     }else if(m==2){
#       surv_curve_data_AS_mod2[[i]] <- surv_curve_data_AS_tmp
#       surv_curve_data_CHQ_mod2[[i]] <- surv_curve_data_CHQ_tmp
#       surv_curve_data_AS_c_mod2[[i]] <- surv_curve_data_AS_c_tmp
#       surv_curve_data_CHQ_c_mod2[[i]] <- surv_curve_data_CHQ_c_tmp
#     }else if(m==3){
#       surv_curve_data_AS_mod3[[i]] <- surv_curve_data_AS_tmp
#       surv_curve_data_CHQ_mod3[[i]] <- surv_curve_data_CHQ_tmp
#       surv_curve_data_AS_c_mod3[[i]] <- surv_curve_data_AS_c_tmp
#       surv_curve_data_CHQ_c_mod3[[i]] <- surv_curve_data_CHQ_c_tmp
#     }
#   }
# }

# save the survival curve data:
# save(surv_curve_data_AS_mod1,surv_curve_data_AS_mod2,surv_curve_data_AS_mod3,surv_curve_data_AS_c_mod1,surv_curve_data_AS_c_mod2,surv_curve_data_AS_c_mod3,
#      surv_curve_data_CHQ_mod1,surv_curve_data_CHQ_mod2,surv_curve_data_CHQ_mod3,surv_curve_data_CHQ_c_mod1,surv_curve_data_CHQ_c_mod2,surv_curve_data_CHQ_c_mod3,
#      file="Simulation-data-mod1-3.RData")

# load the survival curve data:
# load("./Results/Simulation-data-mod1-3.RData") # load survival curve data from the simulations
# Note that the data files for the simulated data were too large to be uploaded to a public repository.
# If you want to reproduce the results, the simulations can be repeated with the provided code or you can send an email to estadler@kirby.unsw.edu.au.

# Plot all survival curves together (Fig. 4C and D, Fig. S12):
m <- 2 # model: 1, 2 or 3
drug <- "CHQ" 
for(q in 1:4){ # for each quartile
  for(i in 1:1000){
    if(m==1){
      if(drug=="AS"){
        data_tmp <- surv_curve_data_AS_mod1[[i]]
        data_tmp_c <- surv_curve_data_AS_c_mod1[[i]]
      }else if(drug=="CHQ"){
        data_tmp <- surv_curve_data_CHQ_mod1[[i]]
        data_tmp_c <- surv_curve_data_CHQ_c_mod1[[i]]
      }
    }else if(m==2){
      if(drug=="AS"){
        data_tmp <- surv_curve_data_AS_mod2[[i]]
        data_tmp_c <- surv_curve_data_AS_c_mod2[[i]]
      }else if(drug=="CHQ"){
        data_tmp <- surv_curve_data_CHQ_mod2[[i]]
        data_tmp_c <- surv_curve_data_CHQ_c_mod2[[i]]
      }
    }else if(m==3){
      if(drug=="AS"){
        data_tmp <- surv_curve_data_AS_mod3[[i]]
        data_tmp_c <- surv_curve_data_AS_c_mod3[[i]]
      }else if(drug=="CHQ"){
        data_tmp <- surv_curve_data_CHQ_mod3[[i]]
        data_tmp_c <- surv_curve_data_CHQ_c_mod3[[i]]
      }
    }
    
    if(q==1){
      if(i==1){
        plot(c(0,data_tmp$time[data_tmp$quartile=='quartile 1'],max(data_tmp$time[data_tmp$quartile=='quartile 1'],data_tmp_c$time[data_tmp_c$quartile=='quartile 1'])),
             c(1,data_tmp$surv[data_tmp$quartile=='quartile 1'],min(data_tmp$surv[data_tmp$quartile=='quartile 1'],data_tmp_c$surv[data_tmp_c$quartile=='quartile 1'])),
             type="s",xlab="Time [days]",ylab = "% Uninfected",las=1,xlim = c(0,400),ylim=c(0,1),col="#7CAF41",main=paste("1000 survival curves for ",
                                                                                                                          drug," treatment (model ",toString(m),")",sep=""))
        points(data_tmp_c$time[data_tmp_c$quartile=='quartile 1'],data_tmp_c$surv[data_tmp_c$quartile=='quartile 1'],pch=3, col="#7CAF41")
      }else{
        lines(c(0,data_tmp$time[data_tmp$quartile=='quartile 1'],max(data_tmp$time[data_tmp$quartile=='quartile 1'],data_tmp_c$time[data_tmp_c$quartile=='quartile 1'])),
              c(1,data_tmp$surv[data_tmp$quartile=='quartile 1'],min(data_tmp$surv[data_tmp$quartile=='quartile 1'],data_tmp_c$surv[data_tmp_c$quartile=='quartile 1'])),
              type="s",col="#7CAF41")
        points(data_tmp_c$time[data_tmp_c$quartile=='quartile 1'],data_tmp_c$surv[data_tmp_c$quartile=='quartile 1'],pch=3, col="#7CAF41")
      }
    }else if(q==2){
      lines(c(0,data_tmp$time[data_tmp$quartile=='quartile 2'],max(data_tmp$time[data_tmp$quartile=='quartile 2'],data_tmp_c$time[data_tmp_c$quartile=='quartile 2'])),
            c(1,data_tmp$surv[data_tmp$quartile=='quartile 2'],min(data_tmp$surv[data_tmp$quartile=='quartile 2'],data_tmp_c$surv[data_tmp_c$quartile=='quartile 2'])),
            type="s",col="#1CBDC2")
      points(data_tmp_c$time[data_tmp_c$quartile=='quartile 2'],data_tmp_c$surv[data_tmp_c$quartile=='quartile 2'],pch=3, col="#1CBDC2")
    }else if(q==3){
      lines(c(0,data_tmp$time[data_tmp$quartile=='quartile 3'],max(data_tmp$time[data_tmp$quartile=='quartile 3'],data_tmp_c$time[data_tmp_c$quartile=='quartile 3'])),
            c(1,data_tmp$surv[data_tmp$quartile=='quartile 3'],min(data_tmp$surv[data_tmp$quartile=='quartile 3'],data_tmp_c$surv[data_tmp_c$quartile=='quartile 3'])),
            type="s",col="#A781BA")
      points(data_tmp_c$time[data_tmp_c$quartile=='quartile 3'],data_tmp_c$surv[data_tmp_c$quartile=='quartile 3'],pch=3, col="#A781BA")
    }else if(q==4){
      lines(c(0,data_tmp$time[data_tmp$quartile=='quartile 4'],max(data_tmp$time[data_tmp$quartile=='quartile 4'],data_tmp_c$time[data_tmp_c$quartile=='quartile 4'])),
            c(1,data_tmp$surv[data_tmp$quartile=='quartile 4'],min(data_tmp$surv[data_tmp$quartile=='quartile 4'],data_tmp_c$surv[data_tmp_c$quartile=='quartile 4'])),
            type="s",col="#F3766E")
      points(data_tmp_c$time[data_tmp_c$quartile=='quartile 4'],data_tmp_c$surv[data_tmp_c$quartile=='quartile 4'],pch=3, col="#F3766E")
    }
  }
}
legend("topright",legend=c("quartile 1", "quartile 2", "quartile 3", "quartile 4"),col=c("#7CAF41","#1CBDC2","#A781BA","#F3766E"), lty = 1, lwd=2)


############################################################################################
### 4. Average number of recurrences per individual and year by time to first recurrence ###
############################################################################################

n_rec_AS_2 <- c() # number recurrences by time to first recurrence for AS and model 2
n_rec_AS_3 <- c() # number recurrences by time to first recurrence for AS and model 3
n_rec_CHQ_2 <- c() # number recurrences by time to first recurrence for CHQ and model 2
n_rec_CHQ_3 <- c() # number recurrences by time to first recurrence for CHQ and model 3

for(l in 1:1000){
  sim2_AS_tmp <- sim2_AS[[l]][[1]]
  sim3_AS_tmp <- sim3_AS[[l]][[1]]
  
  sim2_CHQ_tmp <- sim2_CHQ[[l]][[1]]
  sim3_CHQ_tmp <- sim3_CHQ[[l]][[1]]
  
  for(i in 1:1000){
    if(any(sim2_AS_tmp[,1]==i & sim2_AS_tmp[,2]>0 & sim2_AS_tmp[,4]==1)){# if the ind. has at least one recurrence
      n_tmp <- max(sim2_AS_tmp[sim2_AS_tmp[,1]==i,4])
      n_rec_AS_2 <- rbind(n_rec_AS_2,cbind(sim2_AS_tmp[sim2_AS_tmp[,1]==i & sim2_AS_tmp[,4]==1,2],n_tmp))
    }
    
    if(any(sim3_AS_tmp[,1]==i & sim3_AS_tmp[,2]>0 & sim3_AS_tmp[,4]==1)){# if the ind. has at least one recurrence
      n_tmp <- max(sim3_AS_tmp[sim3_AS_tmp[,1]==i,4])
      n_rec_AS_3 <- rbind(n_rec_AS_3,cbind(sim3_AS_tmp[sim3_AS_tmp[,1]==i & sim3_AS_tmp[,4]==1,2],n_tmp))
    }
    
    if(any(sim2_CHQ_tmp[,1]==i & sim2_CHQ_tmp[,2]>0 & sim2_CHQ_tmp[,4]==1)){# if the ind. has at least one recurrence
      n_tmp <- max(sim2_CHQ_tmp[sim2_CHQ_tmp[,1]==i,4])
      n_rec_CHQ_2 <- rbind(n_rec_CHQ_2,cbind(sim2_CHQ_tmp[sim2_CHQ_tmp[,1]==i & sim2_CHQ_tmp[,4]==1,2],n_tmp))
    }
    
    if(any(sim3_CHQ_tmp[,1]==i & sim3_CHQ_tmp[,2]>0 & sim3_CHQ_tmp[,4]==1)){# if the ind. has at least one recurrence
      n_tmp <- max(sim3_CHQ_tmp[sim3_CHQ_tmp[,1]==i,4])
      n_rec_CHQ_3 <- rbind(n_rec_CHQ_3,cbind(sim3_CHQ_tmp[sim3_CHQ_tmp[,1]==i & sim3_CHQ_tmp[,4]==1,2],n_tmp))
    }
  }
  # print(l)
}

# sort by time to first recurrence:
n_rec_AS_2 <- n_rec_AS_2[order(n_rec_AS_2[,1]),]
n_rec_AS_3 <- n_rec_AS_3[order(n_rec_AS_3[,1]),]
n_rec_CHQ_2 <- n_rec_CHQ_2[order(n_rec_CHQ_2[,1]),]
n_rec_CHQ_3 <- n_rec_CHQ_3[order(n_rec_CHQ_3[,1]),]

# group first recurrence time into intervals of 10 days:
n_rec_AS_2_gr <- matrix(NA,37,1)
n_rec_AS_3_gr <- matrix(NA,37,1)
n_rec_CHQ_2_gr <- matrix(NA,37,1)
n_rec_CHQ_3_gr <- matrix(NA,37,1)

for(i in 1:37){
  tmp <- which(n_rec_AS_2[,1]>(i-1)*10 & n_rec_AS_2[,1]<=i*10)
  n_rec_AS_2_gr[i] <- mean(n_rec_AS_2[tmp,2])
  
  tmp <- which(n_rec_AS_3[,1]>(i-1)*10 & n_rec_AS_3[,1]<=i*10)
  n_rec_AS_3_gr[i] <- mean(n_rec_AS_3[tmp,2])
  
  tmp <- which(n_rec_CHQ_2[,1]>(i-1)*10 & n_rec_CHQ_2[,1]<=i*10)
  n_rec_CHQ_2_gr[i] <- mean(n_rec_CHQ_2[tmp,2])
  
  tmp <- which(n_rec_CHQ_3[,1]>(i-1)*10 & n_rec_CHQ_3[,1]<=i*10)
  n_rec_CHQ_3_gr[i] <- mean(n_rec_CHQ_3[tmp,2])
}

# Visualize (Fig. S13 model estimates): 
plot(seq(5,365,by=10),n_rec_AS_2_gr,type = "l",main="AS, model 2",xlab="Time to first recurrence [days]",ylab="Average number of recurrences",las=1,xlim=c(0,365))
plot(seq(5,365,by=10),n_rec_AS_3_gr,type = "l",main="AS, model 3",xlab="Time to first recurrence [days]",ylab="Average number of recurrences",las=1,xlim=c(0,365))
plot(seq(5,365,by=10),n_rec_CHQ_2_gr,type = "l",main="CHQ, model 2",xlab="Time to first recurrence [days]",ylab="Average number of recurrences",las=1,xlim=c(0,365))
plot(seq(5,365,by=10),n_rec_CHQ_3_gr,type = "l",main="CHQ, model 3",xlab="Time to first recurrence [days]",ylab="Average number of recurrences",las=1,xlim=c(0,365))

