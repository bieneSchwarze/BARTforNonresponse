################################################################################
################################################################################
## LOGIT REGRESSION WITH LASSO REGULARIZATION 
## FOR PREDICTING PARTICIPATION OF TEMPORARY DROPOUT
## Validation on level of individuals (cross validation)
## 07.11.2019
## Sabine Zinn
################################################################################
################################################################################

################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################
library(BART, lib.loc="Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project")
library(glmnet)
rm(list=ls())
set.seed(357)
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

# Add total number of temp. dropout events as (potential) explanatory variable
Dat$numDR <- NA
addNumberTempDR <- function(id){
  st <- Dat$status[Dat$ID_t %in% id]
  st <- cumsum(ifelse(st == 1,1,0))  # only temp. dropout events count
  Dat$numDR[Dat$ID_t %in% id] <<- rep(st[4], 4) # consider only four waves
  return(NULL)
}
nn <- sapply(unique(Dat$ID_t), addNumberTempDR)

# Add sojourn time in state "participation"
sojT <- c()
getSojTime <- function(id){
  #id <- unique(Dat$ID_t)[34]
  vv <- Dat[Dat$ID_t %in% id, "status"]
  s <- rep(NA, length(vv))
  s[1] <- 1
  for(ii in 2:4){
    if(vv[ii]==0){
      s[ii] <- s[ii-1]+1
    } else { 
      s[ii] <- 0
    }
  }
  sojT <<- c(sojT,s)
  return(NULL)
}
nn <- sapply(unique(Dat$ID_t), getSojTime)
Dat$sojT <- sojT
table(Dat$time, Dat$sojT)

################################################################################
# II. DEFINE MODEL MATRICES AND SETTINGS
################################################################################

# Define matrix for times and events (temporary dropout), note that there are no final dropouts in the school context
ids <- unique(Dat$ID_t) 
L <- 4
times <- NULL
delta <- NULL
for(i in 1:length(ids)){
  #i <- 1
  #cat("It: ",i,"\n")
  id <- ids[i]
  mm <- Dat[Dat$ID_t %in% id,]
  if(1 %in% mm$indNV) { 
    fe <- which(mm$indNV %in% 1)[1] 
    mm$indNV[fe:nrow(mm)] <- 1
    wi <- which(mm$indNV %in% 1)
    mm <- mm[-wi,] # drop those already in ind. field & final dropout only occurs in the individual field -> this produces an unbalanced panel    
  }
  if(nrow(mm)>0){ # otherwise: individual moves to individual field in Wave 2 and does not contribute any information to the likelihood on transitions to temp. dropout from Wave 2 onwards
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
table(xx,xy) # there is a slight discrepancy concerning three cases, they seem to be coded wrongly (see below, but okay for the moment)

################################################################################
# II. DEFINE MODEL MATRICES AND SETTINGS
################################################################################
# Add indicator for temporary dropout event
Dat$tempD <- ifelse(Dat$status %in% 1,1,0)

# Define training and test sets (on the level of IDs)
N.train <- round((length(unique(Dat$ID_t))/3)*2)
N.test <- N-N.train
id.train <- sort(sample(x=unique(Dat$ID_t),size=N.train, replace=FALSE)) 
id.test <- sort(setdiff(unique(Dat$ID_t),id.train))

D.train <- Dat[Dat$ID_t %in% id.train,]
trainset <- D.train[,c(91,88,90,89,2:84)]# training set

################################################################################
# III. MODEL ESTIMATION & VALIDATION
################################################################################

# Estimate Logit model (without LASSO)
varE <- apply(trainset, 2,var)  # take out variables without variation 
colnames(trainset)[varE==0] # variable: STATEBE has no variation
trainset <- trainset[,-which(varE==0)]
trainset <- trainset[,-which(colnames(trainset) %in% "STATETH")] # take out "STATETH", model cannot estimate a coefficient for that category
formula <- paste("tempD~",paste(colnames(trainset)[-1], collapse="+"), collapse="")
logit_model <- glm(as.formula(formula), data=trainset, family = binomial(link = "logit"))
summary(logit_model) 

nrow(summary(logit_model)$coefficients)-1 # number coeff is 84 (without intercept)
CI <- confint(logit_model) # warnings occurred "adjusted probability with value zero or one has occurred"
modE <- cbind(coef(logit_model), CI)

# Predict Participation probabilities
D.test <- Dat[!(Dat$ID_t %in% id.train),]
D.test <- D.test[order(D.test$ID_t),]
testset <- D.test[,c(88,90,89,2:84)]
testset <- testset[,!(colnames(testset) %in% c("STATEBE", "STATETH"))] # take out variables without variation 
logit_prob <- predict.glm(logit_model, testset, type = "response")

# Check accuracy, wave-wise
logit_predict <- ifelse(logit_prob>0.5,1,0)
for(wave in 1:4){
  cat("Wave:",wave,"\n")
  print(table(pred=logit_predict[testset$time == wave], obs=D.test$tempD[D.test$time == wave]))
  print(mean(logit_predict[testset$time == wave]==D.test$tempD[D.test$time == wave])) 
}
# Wave 1: 0.94
# Wave 2: 0.99
# Wave 3: 0.99
# Wave 4: 0.98

# Estimate Logit model (with LASSO)
x <- model.matrix(D.train$tempD~., trainset[,-1])
cv.out <- cv.glmnet(x,D.train$tempD, alpha=1, family="binomial", type.measure = "mse") # find optimal smoothin parameter lambda minimizing (total) mean squared error by grid search
plot(cv.out) # optimal lambda (balance accuracy and simplicity: the few coeff as poss but necess)
cv.out$lambda.min # min value of lambda
cv.out$lambda.1se # best vlaue of lambda
coefCV <- coef(cv.out, s=cv.out$lambda.1se) # regression coef: only these coef have non-zero coef (all others have been set to zero)
nrow(summary(coefCV)) # only 8 of 83 coef important
# -> simpler function should be preferred because it is less likely overfitting the data

# Check accuracy, wave-wise
x_test <- model.matrix(D.test$tempD~., testset)
lasso_prob <- predict(cv.out, newx=x_test, s= cv.out$lambda.1se, type="response")
lasso_predict <- ifelse(lasso_prob>0.5,1,0)
for(wave in 1:4){
  cat("Wave:",wave,"\n")
  print(table(pred=lasso_predict[testset$time == wave], obs=D.test$tempD[D.test$time==wave]))
  print(mean(lasso_predict[testset$time == wave]==D.test$tempD[D.test$time==wave])) 
}
# Wave 1: 0.94
# Wave 2: 0.99
# Wave 3: 0.99
# Wave 4: 0.96


