################################################################################
################################################################################
## LOGIT REGRESSION WITH LASSO REGULARIZATION 
## FOR PREDICTING PARTICIPATION OF PERMANENT DROPOUT
## 27.10.2019
## Sabine Zinn
################################################################################
################################################################################


################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################
rm(list=ls())
library(glmnet)
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
xset <- D[, c(88,2:84)]

################################################################################
# III. MODEL ESTIMATION & VALIDATION
################################################################################

# Estimate Logit model (without LASSO)
varE <- apply(xset, 2,var)  # take out variables without variation 
colnames(xset)[varE==0] # variable: STATEBE has no variation
xset <- xset[,-which(varE==0)]
xset <- xset[,-which(colnames(xset) %in% "STATETH")] # take out "STATETH", model cannot estimate a coefficient for that category
formula <- paste("indNV~",paste(colnames(xset), collapse="+"), collapse="")
logit_model <- glm(as.formula(formula), data=D, family = binomial(link = "logit"))
summary(logit_model) # Note: STATETH NA!

nrow(summary(logit_model)$coefficients)-1 # number coeff is 83 (without intercept)
CI <- confint(logit_model)
modE <- cbind(coef(logit_model), CI)

# Predict Participation probabilities
predset <- D[,c(88,2:84)]
predset <- predset[,!(colnames(predset) %in% c("STATEBE", "STATETH"))] # take out variables without variation 
logit_prob <- predict.glm(logit_model, predset, type = "response") # get warnings

# Estimate Logit model (with LASSO)
x <- model.matrix(D$indNV~., xset)
cv.out <- cv.glmnet(x,D$indNV, alpha=1, family="binomial", type.measure = "mse") # find optimal smoothin parameter lambda minimizing (total) mean squared error by grid search
plot(cv.out) # optimal lambda=0.001438362 (balance accuracy and simplicity: the few coeff as poss but necess)
cv.out$lambda.min # min value of lambda
cv.out$lambda.1se # best vlaue of lambda
coefCV <- coef(cv.out, s=cv.out$lambda.1se) # regression coef: only these coef have non-zero coef (all others have been set to zero)
nrow(summary(coefCV)) # only 13 of 83 coef important
# -> simpler function should be preferred because it is less likely overfitting the data

################################################################################
# IV. PLOT IT
################################################################################
load("Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project\\Results\\LASSO\\modelOutput_indF_lasso.RData")
par(mar=c(4, 5, 4, 2), mfrow=c(2,1))
modP <- modE[-1,]
pchS <- ifelse(sign(modP[,2])==sign(modP[,3]), 19, 1)
colS <- c("black", "red", "green", "blue", "orange", "brown")
colP <- c(rep(colS[1],1), rep(colS[2], 27), rep(colS[3],6), 
          rep(colS[4],27), rep(colS[5],7), rep(colS[6],15))
plot(modP[,1], ylab='Effect Size', ylim=c(-3, 7),
     xlab="", pch=pchS, col=colP, axes=FALSE, cex.lab=0.75)
axis(2, cex.axis=0.75)
labs <- rownames(modP)
labs[1] <- "t"
labs <- sub(labs, pattern="TRUE", replacement="")
axis(1, at=1:length(rownames(modP)), labels=labs, las=2, cex.axis=0.75, cex.axis=0.75)
lines(0:82, rep(0,83))
legend("topleft", c( "Wave                                       ",
                     "Student attributes", "Missing value dummies",
                     "Aggregrated student information on school level",
                     "School context information",
                     "Federal state dummies"), fill=colS, ncol=1, 
       cex=0.75, bty="n")

modLASS <- rep(NA, nrow(modE))
indR <- summary(coefCV)[,1]
indR[2:length(indR)] <- (indR[2:length(indR)])-1 # Intercept is double in entry (why ever)
modLASS[indR] <- summary(coefCV)[,3] 
names(modLASS) <- coefCV@Dimnames[[1]][-1]
modL <- cbind(modE, modLASS) 

