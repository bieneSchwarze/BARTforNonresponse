################################################################################
################################################################################
## LOGIT REGRESSION WITH LASSO REGULARIZATION 
## FOR PREDICTING PARTICIPATION OF TEMPORARY DROPOUT
## 07.11.2019
## Sabine Zinn
################################################################################
################################################################################

################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################
library(BART, lib.loc="Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project")
library(glmnet)
library(mice)
rm(list=ls())
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
Dat$tempD <- ifelse(Dat$status %in% 1,1,0)
table(Dat$tempD, Dat$time)
xset <- Dat[,c(91,88,90,89,2:84)]
colnames(xset) <- colnames(Dat[,c(91,88,90,89,2:84)])

################################################################################
# III. MODEL ESTIMATION & VALIDATION
################################################################################

# Estimate Logit model (without LASSO)
varE <- apply(xset, 2,var)  # take out variables without variation 
colnames(xset)[varE==0] # variable: STATEBE has no variation
xset <- xset[,-which(varE==0)]
xset <- xset[,-which(colnames(xset) %in% "STATETH")] # take out "STATETH", model cannot estimate a coefficient for that category
formula <- paste("tempD~",paste(colnames(xset[,-1]), collapse="+"), collapse="")
logit_model <- glm(as.formula(formula), data=xset, family = binomial(link = "logit"))
summary(logit_model) 

nrow(summary(logit_model)$coefficients)-1 # number coeff is 84 (without intercept)
CI <- confint(logit_model) # warnings occurred "adjusted probability with value zero or one has occurred"
modE <- cbind(coef(logit_model), CI)

# Estimate Logit model (with LASSO)
x <- model.matrix(xset$tempD~., xset[,-1])
cv.out <- cv.glmnet(x,xset$tempD, alpha=1, family="binomial", type.measure = "mse") # find optimal smoothin parameter lambda minimizing (total) mean squared error by grid search
plot(cv.out) # optimal lambda (balance accuracy and simplicity: the few coeff as poss but necess)
cv.out$lambda.min # min value of lambda
cv.out$lambda.1se # best vlaue of lambda
coefCV <- coef(cv.out, s=cv.out$lambda.1se) # regression coef: only these coef have non-zero coef (all others have been set to zero)
nrow(summary(coefCV)) # only 11 of 84 coef important
# -> simpler function should be preferred because it is less likely overfitting the data

################################################################################
# IV. PLOT IT
################################################################################
pdf("importExt_complete_lasso_pM.pdf", width=25, height = 10)
par(mar=c(25, 5, 4, 2))
labsN <- c("t","v", "N", "sex", "age", "migration background: yes",
           "number books at home", "native language: German", 
           "percpetual speed", "reasoning", "orthography",
           "reading speed", "mathematical competence", "reading competence",
           "satisfaction with life", "satisfaction with current standards of living",
           "satisfaction with health", "satisfaction with family", "satisfaction with friends",
           "satisfaction with school", "subjective health", "self-esteem",
           "number of people living at home", "recent class repeated: yes",
           "self-concept German", "self-concept mathematics", "self-concept school",
           "grade in German", "grade in mathematics", "number of sick days",
           "self-esteem:NA", "self-concept German:NA", "self-concept mathematics:NA", 
           "self-concept school:NA", "grade in German:NA", "grade in mathematics:NA",
           "proportion females in grade 5", "mean age in in grade 5", 
           "proportion students with migration background in grade 5",  
           "mean number books at home in grade 5", 
           "proportion students with native language German in grade 5",
           "mean percpetual speed in grade 5", "mean reasoning score in grade 5", 
           "mean orthography score in grade 5",
           "mean reading speed in grade 5", "mean mathematical competence in grade 5", 
           "mean reading competence in grade 5",
           "mean satisfaction with life in grade 5", "mean satisfaction with current standards of living in grade 5",
           "mean satisfaction with health in grade 5", "mean satisfaction with family in grade 5", 
           "mean satisfaction with friends in grade 5",
           "mean satisfaction with school in grade 5", "mean subjective health in grade 5", 
           "mean self-esteem in grade 5",
           "mean number of people living at home in grade 5", 
           "proportion of recent class repeaters in grade 5",
           "mean self-concept German in grade 5", "mean self-concept mathematics in grade 5", 
           "mean self-concept in grade 5",
           "mean grade in German in grade 5", "mean grade in mathematics in grade 5", 
           "mean number of sick days in grade 5",
           "school type: Realschule/Gymnasium", "school type: other school", "number of school classes in grade 8 (2013)",
           "number of students in grade 8 (2013)", "funding: public", "city-country: semi-urban", 
           "city-country: rural", "Federal State: Brandenburg", "Federal State: Baden-W.", "Federal State: Bayern",
           "Federal State: Bremen", "Fedaral State: Hessen", "Federal State: Hamburg", "Federal State: Mecklenburg-VP",
           "Fderal State: Niedersachsen", "Federal State: Nordrhein-Westf.", "Federal State: Rheinland-Pfalz",
           "Federal State: Schleswig-Holstein", "Federal State: Saarland", "Federal State: Sachsen", 
           "Federal State: Sachsen-Anhalt", "Federal State: Th?ringen")
modLASS <- rep(0, length(labsN))
modLASS[as.vector(coefCV@i-1)[-1]] <- coefCV@x[-1]
names(modLASS) <- labsN
pchS <- ifelse(modLASS==0,1,19)
#colS <- c("black", "red", "green", "blue", "orange", "brown")
colS <- c("black", gray(0.25), gray(0.35), gray(0.45), gray(0.55), gray(0.65))
colP <- c(rep(colS[1],3), rep(colS[2], 27), rep(colS[3],6), 
          rep(colS[4],27), rep(colS[5],7), rep(colS[6],15))
plot(modLASS, ylab='Effect Size', ylim=c(-6, 8),
     xlab="", pch=pchS, col=colP, axes=FALSE, cex.lab=1.2)
axis(2, cex.axis=0.75)
axis(1, at=1:length(modLASS), labels=labsN, las=2, cex.lab=0.8)
lines(0:85, rep(0,86))
legend("topleft", c("Wave, sojourn time, number previous dropouts                   ",
                    "Student attributes", "Missing value dummies",
                    "Aggregrated student information on school level",
                    "School context information",
                    "Federal state dummies"), fill=colS, ncol=1, 
       cex=1.2, bty="n")
dev.off()

load("Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project\\Results\\partModel_est.Rdata")

pdf("compare_pM.pdf", width=20, height = 10)
par(mar=c(25, 5, 4, 5))
selC <- c(1,2,3,6,7,20,24,30,39,40,45,50,56,65,67,69,71,74,79,80,82)
#cbind(labsN, 1:85)
redC <- post$varprob.mean[selC]
P <-  length(post$varprob.mean)
pchS <- ifelse(post$varprob.mean > 1/P, 19, 1)[selC]
redC[redC < 1/P] <- 0
plot(redC, ylab='Selection Probability', ylim=c(-0.35, 0.35),
     xlab="", pch=pchS, axes=FALSE, cex.lab=1.2, type="p")
axis(2)
labsO <- names(post$varprob.mean)
axis(1, at=1:length(redC), labels=labsN[selC], las=2, cex.lab=0.8)
abline(h=0)
par(new=T)
pchL <- ifelse(modLASS==0,0,15)[selC]
plot(rep(-8,length(selC)),ylab='', ylim=c(-6, 6),
     xlab="", pch=pchL, axes=FALSE, cex.lab=1.2, type="p", col="white")
points(c(1:length(selC))+0.1, modLASS[selC], pch=pchL)
axis(4, ylim=c(-6, 8))
mtext(side=4, line=3, "Effect Size", cex=1.2)
legend("topleft", c("BART: not relevant", "BART: relevant", "LASSO Logit: not relevant", "LASSO Logit: relevant"), 
       pch=c(1,19,0,15), ncol=1, cex=1.2, bty="n")
dev.off()















