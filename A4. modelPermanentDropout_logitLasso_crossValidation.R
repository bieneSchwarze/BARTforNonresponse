################################################################################
################################################################################
## LOGIT REGRESSION WITH LASSO REGULARIZATION 
## FOR PREDICTING PARTICIPATION OF PERMANENT DROPOUT
## Validation on level of individuals (cross-validation)
## 27.10.2019
## Sabine Zinn
################################################################################
################################################################################


################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################
rm(list=ls())
library(glmnet)
set.seed(268)
load("Z:/Projects/p000139_Methoden_Survey_Psych/BART_Project/1. Prepare data/data.Rdata")
N <- length(unique(datl$ID_t))  
time <- rep(1:9, N)
datl <- cbind(datl, time)
Dat <- datl[datl$time %in% 2:5,]
time <- rep(1:4, N)
Dat$time <- time

# # Add total number of participation events as (potential) explanatory variable
# Dat$numPart <- NA
# addNumberParticiaption <- function(id){
#   st <- Dat$status[Dat$ID_t %in% id]
#   st <- cumsum(ifelse(st == 0,1,0))  # only participation events count
#   Dat$numPart[Dat$ID_t %in% id] <<- rep(st[4], 4) # consider only four waves
#   return(NULL)
# }
# nn <- sapply(unique(Dat$ID_t), addNumberParticiaption)

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
trainset <- D.train[,c(88,2:84)] # training set

################################################################################
# III. MODEL ESTIMATION & VALIDATION
################################################################################

# Estimate Logit model (without LASSO)
varE <- apply(trainset, 2,var)  # take out variables without variation 
colnames(trainset)[varE==0] # variable: STATEBE has no variation
trainset <- trainset[,-which(varE==0)]
trainset <- trainset[,-which(colnames(trainset) %in% "STATETH")] # take out "STATETH", model cannot estimate a coefficient for that category
formula <- paste("indNV~",paste(colnames(trainset), collapse="+"), collapse="")
logit_model <- glm(as.formula(formula), data=D.train, family = binomial(link = "logit"))
summary(logit_model) # number coeff
nrow(summary(logit_model)$coefficients)-1 # number coeff is 83 (without intercept)
CI <- confint(logit_model)
modE <- cbind(coef(logit_model), CI)

# Predict Participation probabilities
D.test <- D[!(D$ID_t %in% id.train),]
D.test <- D.test[order(D.test$ID_t),]
testset <- D.test[,c(88,2:84)]
testset <- testset[,!(colnames(testset) %in% c("STATEBE", "STATETH"))] # take out variables without variation 
logit_prob <- predict.glm(logit_model, testset, type = "response")

# Check accuracy, wave-wise
logit_predict <- ifelse(logit_prob>0.5,1,0)
for(wave in 1:4){
  cat("Wave:",wave,"\n")
  print(table(pred=logit_predict[testset$time == wave], obs=D.test$indNV[D.test$time == wave]))
  print(mean(logit_predict[testset$time == wave]==D.test$indNV[D.test$time == wave])) 
}
# Wave 1: 0.94
# Wave 2: 0.89
# Wave 3: 0.91
# Wave 4: 0.90

# Estimate Logit model (with LASSO)
x <- model.matrix(D.train$indNV~., trainset)
cv.out <- cv.glmnet(x,D.train$indNV, alpha=1, family="binomial", type.measure = "mse") # find optimal smoothin parameter lambda minimizing (total) mean squared error by grid search
plot(cv.out) # optimal lambda=0.001438362 (balance accuracy and simplicity: the few coeff as poss but necess)
cv.out$lambda.min # min value of lambda
cv.out$lambda.1se # best vlaue of lambda
coefCV <- coef(cv.out, s=cv.out$lambda.1se) # regression coef: only these coef have non-zero coef (all others have been set to zero)
nrow(summary(coefCV)) # only 8 of 83 coef important
# -> simpler function should be preferred because it is less likely overfitting the data

# Check accuracy, wave-wise
x_test <- model.matrix(D.test$indNV~., testset)
lasso_prob <- predict(cv.out, newx=x_test, s= cv.out$lambda.1se, type="response")
lasso_predict <- ifelse(lasso_prob>0.5,1,0)
for(wave in 1:4){
  cat("Wave:",wave,"\n")
  print(table(pred=lasso_predict[testset$time == wave], obs=D.test$indNV[D.test$time==wave]))
  print(mean(lasso_predict[testset$time == wave]==D.test$indNV[D.test$time==wave])) 
}
# Wave 1: 0.94
# Wave 2: 0.89
# Wave 3: 0.91
# Wave 4: 0.89
