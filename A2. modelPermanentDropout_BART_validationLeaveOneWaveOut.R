################################################################################
################################################################################
## BART FOR PREDICTING PERMANENT DROPOUT
## Validation on level of waves (validation 'leave last wave out')
## 05.11.2019
## Sabine Zinn
################################################################################
################################################################################

################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################
rm(list=ls())
set.seed(2567)
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
# 1. Define training and test sets (on the level of waves)
D.train <- D[D$time %in% 1:3,]
x.train  <- D.train[!duplicated(D.train$ID_t), c(2:84)]
colnames(x.train) <- colnames(D[,c(2:84)]) 

L <- 3
times <- rep(NA, N)
delta <- rep(NA, N)
ids <- unique(D.train$ID_t)
for(i in 1:length(ids)){
  id <- ids[i]
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
save.image("Z:/Projects/p000139_Methoden_Survey_Psych/BART_Project/Results/modelIndF_est_validLevelWave.Rdata")

################################################################################
# V. VALIDATION on the level of wave 4
# (in paper this refers to wave 5, here counting starts at wave 0)
################################################################################
D.test <- D[D$time %in% c(3,4),] 
  # take in addition wave 3, otherwise the prediction function of the BART package does not work;
  # however, only consider predictions for wave 4 (this is just a technical trick)
testset  <- D.test[, c(1,88,2:84, 87)]
colnames(testset) <- colnames(D[,c(1,88,2:84,87)]) 
x.test <- testset[,!(colnames(trainset) %in% c("ID_t", "indNV",colnames(post$tx.train)[-(post$rm.const)]))]
colnames(x.test)[1] <- "t"
predV <- predict(post, x.test) 
head(predV$prob.test[,1:20])
x.test_W4 <- x.test[x.test$t %in% 4,]
predV_W4 <- predV$prob.test[,x.test$t %in% 4]# Take only those predicted values that refer to Wave 4

# Derive posterior draws of prob. of perm. dropout for test set, for wave 4
dr4 <- apply(predV_W4,2,mean) # Wave 4, predicted probability to dropout permanently
bart_predict_wave4 <- ifelse(dr4>0.5,1,0) # accuracy wave 4
print(mean(bart_predict_wave4==D.test$indNV[D.test$time==4])) # 0.90
