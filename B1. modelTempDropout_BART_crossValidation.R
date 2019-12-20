################################################################################
################################################################################
## BART FOR PREDICTING PROBABILITIES OF TEMPORARY DROPOUT
## Validation on level of individuals (cross-validation)
## 06.11.2019
## Sabine Zinn
################################################################################
################################################################################

rm(list=ls())
set.seed(235)
################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################
library(BART, lib.loc="Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project")
library(foreign)
library(mice)
setwd("Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project\\1. Prepare data")
load("data.RData")
N <- length(unique(datl$ID_t)) # N=4559 
time <- rep(1:9, N)
datl <- cbind(datl, time)
Dat <- datl[datl$time %in% 2:5,]
time <- rep(1:4, N)
Dat$time <- time
id.risk <- Dat[Dat$wave==2 & Dat$indNV==0 ,1] # only take those students into the risk set who did not perm. dropout already in wave 2
Dat <- Dat[Dat$ID_t %in% id.risk,]
N <- length(unique(Dat$ID_t)) # N=4266
md.pattern(Dat, plot=FALSE) # very few missing values in books & mig

################################################################################
# II. DEFINE MODEL MATRICES AND SETTINGS
################################################################################
# Matrix of covariates as training data
# 1. Define training and test sets (on the level of IDs)
N.train <- round((N/3)*2) # N=2844
N.test <- N-N.train # N=1422
id.train <- sort(sample(x=unique(Dat$ID_t),size=N.train, replace=FALSE)) 
id.test <- sort(setdiff(unique(Dat$ID_t),id.train))

# Define matrix for times and events, i.e. temporary dropout
L <- 4
times <- NULL
delta <- NULL
for(i in 1:length(id.train)){
  #i <- 1
  #cat("It: ",i,"\n")
  id <- id.train[i]
  mm <- Dat[Dat$ID_t %in% id,]
  if(1 %in% mm$indNV) { # drop those already dopped out permanently -> this produces an unbalanced panel 
     fe <- which(mm$indNV %in% 1)[1] 
     mm$indNV[fe:nrow(mm)] <- 1
     wi <- which(mm$indNV %in% 1)
     mm <- mm[-wi,]    
   }
  if(nrow(mm)>0){ # count events 
    time.i <- 1:nrow(mm)
    delta.i <- rep(0, nrow(mm))
    if(1 %in% mm$status){ # there is an event / are events "temp. dropout"
     td.i <- which(mm$status %in% 1)  
     delta.i[td.i] <- 1
    }
    if(length(time.i)<L)
      time.i <- c(time.i, rep(0, L-length(time.i)))
    if(length(delta.i)<L)
      delta.i <- c(delta.i, rep(0, L-length(delta.i)))
    times <- rbind(times, c(id, time.i))
    delta <- rbind(delta, c(id, delta.i))
    rm(time.i); rm(delta.i)
  } 
}
apply(delta[,-1],2,sum)
xx <- apply(times[,-1],1,paste, collapse="-")
xy <- apply(delta[,-1],1,paste, collapse="-")
table(xx,xy) # empirical distribution of events is okay
sum(apply(delta[,-1],2,sum)) # this is the number of temporary dropouts in school in the traning data set (N=503)

# Define training data set
x.train <- Dat[Dat$ID_t %in% id.train,]
x.train <- x.train[!duplicated(x.train$ID_t), c(2:84)]
colnames(x.train) <- colnames(Dat[,c(2:84)])

# Define init structure for BART
pre <- recur.pre.bart(times=times[,-1], delta=delta[,-1], x.train= as.matrix(x.train))

# Fit model
post <- recur.bart(y.train=pre$y.train, pre$tx.train, x.test=pre$tx.test, 
                    #  ntree = 300, usequants=TRUE, nskip=200, sparse=TRUE, 
                    #  ndpost=1500, keepevery=200)
                    ntree = 300, usequants=TRUE, nskip=500, sparse=TRUE, 
                    ndpost= 5000, keepevery=500)
save.image("Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project\\Results\\partModel_est_validID.Rdata")

################################################################################
# V. VALIDATION (on the level of ids)
################################################################################
# Define matrix for times and events, i.e. temporary dropout
L <- 4
times <- NULL
delta <- NULL
for(i in 1:length(id.test)){
  #i <- 1
  #cat("It: ",i,"\n")
  id <- id.test[i]
  mm <- Dat[Dat$ID_t %in% id,]
  if(1 %in% mm$indNV) { # drop those already dopped out permanently -> this produces an unbalanced panel 
    fe <- which(mm$indNV %in% 1)[1] 
    mm$indNV[fe:nrow(mm)] <- 1
    wi <- which(mm$indNV %in% 1)
    mm <- mm[-wi,]    
  }
  if(nrow(mm)>0){ # count events 
    time.i <- 1:nrow(mm)
    delta.i <- rep(0, nrow(mm))
    if(1 %in% mm$status){ # there is an event / are events "temp. dropout"
      td.i <- which(mm$status %in% 1)  
      delta.i[td.i] <- 1
    }
    if(length(time.i)<L)
      time.i <- c(time.i, rep(0, L-length(time.i)))
    if(length(delta.i)<L)
      delta.i <- c(delta.i, rep(0, L-length(delta.i)))
    times <- rbind(times, c(id, time.i))
    delta <- rbind(delta, c(id, delta.i))
    rm(time.i); rm(delta.i)
  } 
}
apply(delta[,-1],2,sum)
xx <- apply(times[,-1],1,paste, collapse="-")
xy <- apply(delta[,-1],1,paste, collapse="-")
table(xx,xy) # empirical distribution of events is okay
sum(apply(delta[,-1],2,sum))

# Define training data set
x.test <- Dat[Dat$ID_t %in% id.test,]
x.test <- x.test[!duplicated(x.test$ID_t), c(2:84,86)] # CHECK THIS!
colnames(x.test) <- colnames(Dat[,c(2:84,86)])

# Define init structure for BART
pre <- recur.pre.bart(times=times[,-1], delta=delta[,-1], x.train = as.matrix(x.test))
preX <- pre$tx.train[,post$rm.const]
preXP <- preX[,!(colnames(preX) %in% "status")] 
predV <- predict(post,preXP)   # on the observed
head(predV$prob.test[,1:20])

# Derive posterior draws of prob. of temp. dropout for test set, for each wave separately
pi.1V <- predV$prob.test[,preX[,"t"]==1]
pi.2V <- predV$prob.test[,preX[,"t"]==2] 
pi.3V <- predV$prob.test[,preX[,"t"]==3] 
pi.4V <- predV$prob.test[,preX[,"t"]==4] 

# Get accuracy for each wave
preX <- as.data.frame(preX)
preX$tempD <- ifelse(preX$status %in% 1,1,0)
prop.table(table(preX[,"tempD"], preX[,"t"]),2) # empirical (observed) distribution of temporarily dropout events
table(preX[,"t"]) # risk set decreases from wave to wave
dr1 <- mean(apply(pi.1V,1,mean)) # Wave 1, predicted probability to dropout temporarily
bart_predict_wave1 <- ifelse(apply(pi.1V,2,mean)>0.5,1,0) # accuracy wave 1
print(mean(bart_predict_wave1==preX[,"tempD"][preX[,"t"]==1])) # 0.95
dr2 <- mean(apply(pi.2V,1,mean)) # Wave 2, predicted probability to dropout temporarily
bart_predict_wave2 <- ifelse(apply(pi.2V,2,mean)>0.5,1,0) # accuracy wave 2
print(mean(bart_predict_wave2==preX[,"tempD"][preX[,"t"]==2])) # 0.96
dr3 <- mean(apply(pi.3V,1,mean)) # Wave 3, predicted probability to ropout temporarily
bart_predict_wave3 <- ifelse(apply(pi.3V,2,mean)>0.5,1,0) # accuracy wave 3
print(mean(bart_predict_wave3==preX[,"tempD"][preX[,"t"]==3])) # 0.96
dr4 <- mean(apply(pi.4V,1,mean)) # Wave 4, predicted probability to dropout temporarily
bart_predict_wave4 <- ifelse(apply(pi.4V,2,mean)>0.5,1,0) # accuracy wave 4
print(mean(bart_predict_wave4==preX[,"tempD"][preX[,"t"]==4])) # 0.96
