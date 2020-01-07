################################################################################
################################################################################
## BART FOR PREDICTING PERMANENT DROPOUT
## Validation on level of individuals (i.e., cross-validation)
## 28.10.2019
## Sabine Zinn
################################################################################
################################################################################

################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################
rm(list=ls())
set.seed(986)
library(BART, lib.loc="Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project")
load("Z:/Projects/p000139_Methoden_Survey_Psych/BART_Project/1. Prepare data/data.Rdata")
N <- length(unique(datl$ID_t))  
time <- rep(1:9, N)
datl <- cbind(datl, time)
Dat <- datl[datl$time %in% 2:5,]
time <- rep(1:4, N)
Dat$time <- time

# Create unbalanced panel for transition to individual field
table(Dat$indNV, Dat$status) #fin. dropouts only in individual field -> have to model only those
ids <- unique(Dat$ID_t)
D <- NULL
for(i in 1:length(ids)){
 if(i %in% seq(from=1, to=length(ids), by=200))
   cat("It: ",i,"\n")
 m <- Dat[Dat$ID_t %in% ids[i],]
 if(1 %in% m$indNV){
   m <- m[1:which(m$indNV==1)[1],]
 }
 D <- rbind(D,m)
}
table(D$status, D$indNV)
table(D$time, D$indNV)

################################################################################
# II. DEFINE MODEL MATRICES AND SETTINGS
################################################################################
# 1. Define training and test sets (on the level of IDs)
N.train <- round((N/3)*2)
N.test <- N-N.train
id.train <- sort(sample(x=unique(Dat$ID_t),size=N.train, replace=FALSE)) 
id.test <- sort(setdiff(unique(Dat$ID_t),id.train))
  
D.train <- D[D$ID_t %in% id.train,]
x.train <- D.train[!duplicated(D.train$ID_t), c(1:84)]
x.train <- x.train[x.train$ID_t %in% id.train,-1]
colnames(x.train) <- colnames(D[,c(2:84)]) 

L <- 4
times <- rep(NA, N.train)
delta <- rep(NA, N.train)
for(i in 1:N.train){
  id <- id.train[i]
  mm <- D.train[D.train$ID_t %in% id,]
  if( 1 %in% mm$indNV){
    times[i] <- mm$time[which(mm$indNV==1)]
    delta[i] <- 1
  } else{
    #cat("It: ",i,"\n")
    times[i] <- mm$time[nrow(mm)]
    delta[i] <- 0
  }
}
table(times,delta, exclude=NULL)

################################################################################
# III. MODEL ESTIMATION
################################################################################
post <- surv.bart(x.train, times=times, delta=delta,
                  ntree = 300, usequants=TRUE, nskip=500, sparse=TRUE, 
                  ndpost=5000, keepevery=500)
save.image("Z:/Projects/p000139_Methoden_Survey_Psych/BART_Project/Results/modelIndF_est_validLevelIDs.Rdata")

################################################################################
# V. VALIDATION (on the level of ids)
################################################################################
D.test <- D[D$ID_t %in% id.test,]
x.test <- D.test[, c(88,2:84)] # in long format, without ID but with time as predictor
colnames(x.test) <- colnames(post$tx.train) 
predV <- predict(post, x.test[,post$rm.const])
head(predV$prob.test[,1:20])

# Derive posterior draws of prob. of perm. dropout for test set, for each wave separately
pi.1V <- predV$prob.test[,x.test[,"t"]==1]
pi.2V <- predV$prob.test[,x.test[,"t"]==2] 
pi.3V <- predV$prob.test[,x.test[,"t"]==3] 
pi.4V <- predV$prob.test[,x.test[,"t"]==4] 

# Get accuracy for each wave 
prop.table(table(D.test$indNV, D.test$time),2) # empirical (observed) distribution of permanent dropout events
apply(table(D.test$indNV, D.test$time),2,sum) # risk set decreases from wave to wave
dr1 <- mean(apply(pi.1V,1,mean)) # Wave 1, predicted probability to dropout permanently 
bart_predict_wave1 <- ifelse(apply(pi.1V,2,mean)>0.5,1,0) # accuracy wave 1
print(mean(bart_predict_wave1==D.test$indNV[D.test$time==1])) # 0.95
dr2 <- mean(apply(pi.2V,1,mean)) # Wave 2, predicted probability to dropout permanently
bart_predict_wave2 <- ifelse(apply(pi.2V,2,mean)>0.5,1,0) # accuracy wave 2
print(mean(bart_predict_wave2==D.test$indNV[D.test$time==2])) # 0.91
dr3 <- mean(apply(pi.3V,1,mean)) # Wave 3, predicted probability to dropout permanently
bart_predict_wave3 <- ifelse(apply(pi.3V,2,mean)>0.5,1,0) # accuracy wave 3
print(mean(bart_predict_wave3==D.test$indNV[D.test$time==3])) # 0.93
dr4 <- mean(apply(pi.4V,1,mean)) # Wave 4, predicted probability to dropout permanently
bart_predict_wave4 <- ifelse(apply(pi.4V,2,mean)>0.5,1,0) # accuracy wave 4
print(mean(bart_predict_wave4==D.test$indNV[D.test$time==4])) # 0.90
