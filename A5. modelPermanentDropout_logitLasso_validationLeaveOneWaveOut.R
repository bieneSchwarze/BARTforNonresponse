################################################################################
################################################################################
## LOGIT REGRESSION WITH LASSO REGULARIZATION 
## FOR PREDICTING PARTICIPATION OF PERMANENT DROPOUT
## Validation on level of waves (validation 'leave last wave out')
## 28.10.2019
## Sabine Zinn
################################################################################
################################################################################

################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################
rm(list=ls())
library(glmnet)
set.seed(286)
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
# Take last wave as test set, data on remaining waves constitute the training set.
D.train <- D[D$time %in% 1:3,]
trainset  <- D.train[, c(88,2:84, 87)]
colnames(trainset) <- colnames(D[,c(88,2:84, 87)]) 

################################################################################
# III. MODEL ESTIMATION & VALIDATION
################################################################################

# Estimate Logit model (without LASSO)
varE <- apply(trainset, 2,var)  # take out variables without variation 
colnames(trainset)[varE==0] # variable: STATEBE has no variation
trainset <- trainset[,-which(varE==0)]
formula <- paste("indNV~",paste(colnames(trainset)[-84], collapse="+"), collapse="") # CHECK -85
logit_model <- glm(as.formula(formula), data=D.train, family = binomial(link = "logit"))
summary(logit_model) # number coeff
trainset <- trainset[,-which(colnames(trainset) %in% "STATETH")] # take out "STATETH", not enough variation
formula <- paste("indNV~",paste(colnames(trainset)[-83], collapse="+"), collapse="") # CHECK -84
logit_model <- glm(as.formula(formula), data=D.train, family = binomial(link = "logit"))
summary(logit_model) # number coeff
nrow(summary(logit_model)$coefficients)-1 # number coeff is 83 (without intercept)
CI <- confint(logit_model)
modE <- cbind(coef(logit_model), CI)

# Predict Participation probabilities
D.test <- D[D$time==4,]
testset <- D.test[,c(88,2:84)]
testset <- testset[,!(colnames(testset) %in% c("STATEBE", "STATETH"))] # take out variables without variation 
logit_prob <- predict.glm(logit_model, testset, type = "response")

# Check accuracy
logit_predict <- ifelse(logit_prob>0.5,1,0)
table(pred=logit_predict, obs=D.test$indNV[D.test$time == 4])
mean(logit_predict==D.test$indNV[D.test$time == 4]) # 0.90

# Estimate Logit model (with LASSO)
x <- model.matrix(D.train$indNV~., trainset[,-83]) 
cv.out <- cv.glmnet(x,D.train$indNV, alpha=1, family="binomial", type.measure = "mse") # find optimal smoothin parameter lambda minimizing (total) mean squared error by grid search
plot(cv.out) # balance accuracy and simplicity: the few coeff as poss but necess
cv.out$lambda.min # min value of lambda
cv.out$lambda.1se # best vlaue of lambda
coefCV <- coef(cv.out, s=cv.out$lambda.1se) # regression coef: only these coef have non-zero coef (all others have been set to zero)
nrow(summary(coefCV)) # only 10 of 83 coef important (84 if time is seen as a predictor)
# -> simpler function should be preferred because it is less likely overfitting the data

# Check accuracy
x_test <- model.matrix(D.test$indNV~., testset)
lasso_prob <- predict(cv.out, newx=x_test, s= cv.out$lambda.1se, type="response")
lasso_predict <- ifelse(lasso_prob>0.5,1,0)
table(pred=lasso_predict, obs=D.test$indNV[D.test$time==4])
mean(lasso_predict==D.test$indNV[D.test$time==4]) # 0.90



