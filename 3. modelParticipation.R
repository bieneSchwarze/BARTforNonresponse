################################################################################
################################################################################
## BART FOR PREDICTING PROBABILITIES OF TEMPORARY DROPOUT
## Level: school context
## 19.06.2019
## Sabine
################################################################################
################################################################################


################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################
library(BART)
library(foreign)
library(mice)
rm(list=ls())
load("data.RData") # load data prepared using "1. loadData.R"
datL <- datl
datL <- datL[datL$wave %in% 1:5,]
datL <- datL[order(datL$ID_t),]    # very few missing values in books & mig
md.pattern(datL, plot=FALSE)

################################################################################
# II. DEFINE MODEL MATRICES AND SETTINGS
################################################################################
# Matrix of covariates as training data
x.train <- datL[!duplicated(datL$ID_t),2:84]
colnames(x.train) <- colnames(datL[,2:84])
# Define matrix for times and events (temporary dropout), note that there are no final dropouts in the school context
ids <- unique(datL$ID_t) # N=4559
L <- 4
times <- NULL
delta <- NULL
for(i in 1:length(ids)){
  id <- ids[i]
  mm <- datL[datL$ID_t %in% id,][-1,]
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
table(xx,xy) 
id2 <- datL[datL$wave==2 & datL$indNV==0,1] # use only those cases that are at risk in wave 2

# Define init structure for BART
pre <- recur.pre.bart(times=times[,-1], delta=delta[,-1], x.train= as.matrix(x.train[ids %in% id2,]))

# Define test data to derive Friedman's partial dependence for single varaibles 
NK <- nrow(pre$tx.test)
pre$tx.test <- rbind(pre$tx.test, pre$tx.test, pre$tx.test, pre$tx.test)
pre$tx.test[1:NK,"sex"] <- rep(0, NK)             ## sex=0
pre$tx.test[(NK+1):(2*NK), "sex"] <- rep(1, NK)   ## sex=1
pre$tx.test[(2*NK+1):(3*NK), "mig"] <- rep(0, NK) ## mig=0
pre$tx.test[(3*NK+1):(4*NK), "mig"] <- rep(1, NK) ## mig=1

# Fit model
set.seed(135)
post <- recur.bart(y.train=pre$y.train, pre$tx.train, x.test=pre$tx.test, 
                    ntree = 300, usequants=TRUE, nskip=500, sparse=TRUE, 
                    ndpost= 5000, keepevery=500)

################################################################################
# III. CHECK CONVERGENCE
################################################################################
# A. There is some autocorrelation: 
# (1) Close chain elements are very similar, meaning that throwing one away does not loose a lot of information (that is what the autocorrelation plot shows)
# (2) You need a lot of repetitions to get convergence, meaning that you get very large chains if you don't thin. Because of that, working with the full chain can be very slow, 
#     costs a lot of storage, or even lead to memory problems when monitoring a lot of variables.
N <- ncol(post$yhat.train) / pre$K 
k <- floor(seq(1, N, length.out=5))
j. <- seq(-0.5, 0.4, length.out=5)
K <- pre$K
par(mfrow = c(2,2))
for(j in 1:K) {
  for(i in c(1:5)) {
    h <- (k[i]-1)*K+j
    auto.corr <- acf(post$yhat.test[ , h], plot=FALSE)
    max.lag <- max(auto.corr$lag[ , 1, 1])
    if(i==1)
      plot(1:max.lag+j.[i], auto.corr$acf[1+(1:max.lag), 1, 1],
           type='h', xlim=c(0, max.lag+1), ylim=c(-1, 1),
           sub=paste0('t=', pre$times[j]), ylab='Autokorrelation', xlab='lag')
    else
      lines(1:max.lag+j.[i], auto.corr$acf[1+(1:max.lag), 1, 1],
            type='h', col=i)
  }
  Sys.sleep(1)
}

# B. Trace plot: traces demonstrate that samples of f(x_{ij}) do not adaequately traverse sample space. 
par(mfrow = c(2,2), mar=c(6, 5, 2, 3))
k <- floor(seq(1, N, length.out=10))
for(j in 1:K) {
  for(i in 1:10) {
    h <- (k[i]-1)*K+j
    if(i==1)
      plot(post$yhat.test[ , h], type='l',
           ylim=c(-4, 4), sub=paste0('t=', pre$times[j]),
           ylab=expression(Phi(g(x, T, M))), xlab='Züge')
    else
      lines(post$yhat.test[ , h], type='l', col=i)
  }
  Sys.sleep(1)
}

# C. Geweke convergence plot: diagnostic for probit BART: 
LL <- 10
k <- floor(seq(1, N, length.out=LL))
K <- pre$K
geweke <- list(1:LL)
for(i in 1:LL) {
  h <- (k[i]-1)*K+1:K
  geweke[[i]] <- gewekediag(post$yhat.test[ ,h])
}
max.t <- max(pre$times)
min.t <- -max.t/LL
for(i in 1:LL) {
  if(i==1) {
    plot(pre$times, geweke[[i]]$z, type='p',
         ylab='z-Score', xlab='t', ylim=c(-5, 5), xlim=c(min.t, max.t))
    lines(pre$times, rep(-1.96, K), type='l', col=6)
    lines(pre$times, rep(+1.96, K), type='l', col=6)
    lines(pre$times, rep(-2.576, K), type='l', col=5)
    lines(pre$times, rep(+2.576, K), type='l', col=5)
    lines(pre$times, rep(-3.291, K), type='l', col=4)
    lines(pre$times, rep(+3.291, K), type='l', col=4)
    lines(pre$times, rep(-3.891, K), type='l', col=3)
    lines(pre$times, rep(+3.891, K), type='l', col=3)
    lines(pre$times, rep(-4.417, K), type='l', col=2)
    lines(pre$times, rep(+4.417, K), type='l', col=2)
    text(c(0, 0), c(-1.96, 1.96), pos=2, cex=0.6, labels='0.95')
    text(c(0, 0), c(-2.576, 2.576), pos=2, cex=0.6, labels='0.99')
    text(c(0, 0), c(-3.291, 3.291), pos=2, cex=0.6, labels='0.999')
    text(c(0, 0), c(-3.891, 3.891), pos=2, cex=0.6, labels='0.9999')
    text(c(0, 0), c(-4.417, 4.417), pos=2, cex=0.6, labels='0.99999')
  }
  else points(pre$times, geweke[[i]]$z)#type='l')
}

################################################################################
# IV. VARIABLE SELECTION & TREES 
################################################################################
################################################################################
# A. Variable importance
par(mar=c(7, 5, 4, 2))
P <-  length(post$varprob.mean)
pchS <- ifelse(post$varprob.mean > 1/P, 19, 1)
colS <- c("black", "red", "green", "blue", "orange", "brown")
colP <- c(rep(colS[1],3), rep(colS[2], 27), rep(colS[3],6), 
          rep(colS[4],27), rep(colS[5],7), rep(colS[6],15))
plot(post$varprob.mean, ylab='Selection Probability', ylim=c(0, 0.6),
     xlab="", pch=pchS, col=colP, axes=FALSE)
axis(2)
axis(1, at=1:P, labels=names(post$varprob.mean), las=2)
lines(c(0, length(names(post$varprob.mean))), c(1/P, 1/P))
legend("topright", c("Wave, sojourn time, number previous dropouts",
                     "Student attributes", "Missing value dummies",
                     "Aggregrated student information on school level",
                     "School context information",
                     "Federal state dummies"), fill=colS, ncol = 1,
                      cex=0.9, bty="n")

# B. Trees
tc <- textConnection(post$treedraws$tree)
trees <- read.table(file=tc, fill=TRUE,
                    row.names=NULL, header=FALSE,
                    col.names=c('node', 'var', 'cut', 'leaf'))
close(tc)
m <- 1 ## MCMC sample, count it
M <- trees$node[1] # number of max. iterations 
n <- 0 ## trees
H <- trees$var[1] # 200 trees
branch <- matrix(0, nrow=P, ncol=P)
dimnames(branch)[[1]] <- paste0('x', 1:P)
dimnames(branch)[[2]] <- paste0('x', 1:P)
L <- nrow(trees)
for(l in 2:L) {
    if(is.na(trees$leaf[l])) {
        n <- n+1
        if(n>H) {
            n <- 1
            m <- m+1
        }
        C <- trees$node[l] ## nodes in tree
        B <- (C-1)/2 ## branches in tree
        i <- 0
        j <- 0
        if(C>1) vars <- integer(C)
        branch. <- 0*branch
    }
    else if(B>1) {
        i <- i+1
        h <- trees$node[l]
        if(i<C) {
            t <- floor(log2(h))
            k <- h-2^t
            if(trees$node[l+1]==(2^(t+1)+2*k)) {
                vars[h] <- trees$var[l]+1
                j <- j+1
                if(j>B) stop('Too many branches')
            }
        }
        else {
            for(h. in (C-1):2) {
                h <- h.
                j <- vars[h]
                if(j!=0) {
                    for(t in (floor(log2(h))-1):0) {
                        if((h%%2)==0) k <- (h-2^(t+1))/2
                        else k <- (h-2^(t+1)-1)/2
                        h <- 2^t+k
                        i <- vars[h]
                        if(i!=j) branch.[min(i, j), max(i, j)] <- 1
                        vars[h] <- 0
                    }
                }
            }
            branch <- branch+branch.
        }
    }
}
C <- sum(c(branch))
for(i in 1:(P-1))
    for(j in (i+1):P)
        if(i!=j) branch[j, i] <- branch[i, j]/C
XX <- round(branch, digits=2)
colnames(XX)=rownames(XX)=names(post$varprob.mean)
xx <- as.vector(XX)
range(xx[xx<1])
ii <- which(XX>=0.02 & XX<1, arr.ind=TRUE)
RES <- cbind(colnames(XX)[ii[,2]], rownames(XX)[ii[,1]], XX[ii]) # most frequent splits
RES[order(RES[,3], RES[,1], decreasing = TRUE),]
colnames(RES) <- c("firstNode", "followingNode", "relFreqOfallBranches")

table(is.na(trees[-1,2])) # number of trees in total (num interations x times of trees per iteration)
tab <- table(trees[is.na(trees[,2]), 1]) # distribution of total number of nodes in all trees 
round(tab / sum(tab),3)

tab <- table(trees[which(is.na(trees[,2]))+1,2]) # freq. distr. of first node variable
names(tab) <- names(post$varcount.mean)
tab <- sort(tab, decreasing=TRUE)
round(tab / sum(tab),3)

#head(trees) # have a look at one tree
#names(post$varcount.mean)[81]
#post$treedraws$cutpoints$STATESH[1]
#post$treedraws$cutpoints$t[1]

################################################################################
# V. PREDICTION
################################################################################
ppO <- predict(post, pre$tx.train[,post$rm.const])   # on the observed

# A. Relative risk via Friedman's partial dependence function
M <- nrow(post$yhat.test)
RI.male.female <- matrix(0, nrow=M, ncol=pre$K)
RI.nmig.mig <- matrix(0, nrow=M, ncol=pre$K)
for(j in 1:K) {
    h <- seq(j, NK, K)
    RI.male.female[ , j] <- apply(post$prob.test[,1:NK][,h]/
                            post$prob.test[,(NK+1):(2*NK)][,h], 1, mean)
    RI.nmig.mig[ , j] <- apply(post$prob.test[,(2*NK+1):(3*NK)][,h]/
                            post$prob.test[,(3*NK+1):(4*NK)][,h], 1, mean)
}
RI.male.female.mu <- apply(RI.male.female, 2, mean)
RI.male.female.025 <- apply(RI.male.female, 2, quantile, probs=0.025)
RI.male.female.975 <- apply(RI.male.female, 2, quantile, probs=0.975)
RI.nmig.mig.mu <- apply(RI.nmig.mig, 2, mean)
RI.nmig.mig.025 <- apply(RI.nmig.mig, 2, quantile, probs=0.025)
RI.nmig.mig.975 <- apply(RI.nmig.mig, 2, quantile, probs=0.975)
par(mfrow=c(2,1))
plot(post$times, RI.male.female.mu, col='blue',
     log='y', main='Participation Propensity: Male vs. Female',
     type='l', ylim=c(0.9, 1.1), ylab='RI(t)', xlab="",
     axes=FALSE)
axis(1, at=1:4, labels=c("Wave 2", "Wave 3", "Wave 4", "Wave 5"))
axis(2)
lines(post$times, RI.male.female.025, lty=2, col="blue")
lines(post$times, RI.male.female.975, lty=2, col="blue")
plot(post$times, RI.nmig.mig.mu, col='blue',
     log='y', main='Participation Propensity: NonMig vs. Mig',
     type='l', ylim=c(0.9, 1.1), ylab='RI(t)', xlab="",
     axes=FALSE)
axis(1, at=1:4, labels=c("Wave 2", "Wave 3", "Wave 4", "Wave 5"))
axis(2)
lines(post$times, RI.nmig.mig.025, lty=2, col='blue')
lines(post$times, RI.nmig.mig.975, lty=2, col='blue')

# B. Predict probabilities of temp. dropout
p.1 <- ppO$prob.test[,pre$tx.train[,"t"] %in% 1] 
hist(p.1)
median(p.1); mean(p.1)
p.2 <- ppO$prob.test[,pre$tx.train[,"t"] %in% 2]
hist(p.2)
median(p.2); mean(p.2)
p.3 <- ppO$prob.test[,pre$tx.train[,"t"] %in% 3]
hist(p.3)
median(p.3); mean(p.3)
p.4 <- ppO$prob.test[,pre$tx.train[,"t"] %in% 4]
hist(p.4)
median(p.4); mean(p.4)
tdList <- list(p1=p.1, p2=p.2, p3=p.3, p4=p.4)

# Function to derive the highest posterior density
getHPD <- function(thetas, alpha=0.1){
  thetas <- sort(thetas)
  n <- length(thetas)
  bo <- trunc((1-alpha)*n)
  int_w <- c(thetas[1],thetas[1+bo])
  di <- diff(int_w)  
  sm <- 1
  for(j in 1:(n-bo)){
    #cat("j", j, "\n")
    int_j <- c(thetas[j],thetas[j+bo])
    diff_j <- diff(int_j)
    if(diff_j < di) {
      di <- diff_j
      sm <- j
      int_w <- int_j
    }
  }
  return(list(int_w=int_w, wj=sm))
}

# Have a look at posterior draws for prob. to temp. dropout for ten indiv. at wave 2
par(mfrow=c(4,3))
piL <- tdList[[1]]
for(ind in 1:10){ # for all individuals (here n=10) -> one distribution for the estimated prob (i.e. one distribution per ind.)
  dd <- density(piL[,ind])
  plot(x=dd$x, y=dd$y, lty=1, type="l", col="grey35")
  points(median(piL[,ind]),0)
}

# Derive prob. of temp. dropout, averaged over all individuals
emp <- 1-c(0.953, 0.947, 0.946, 0.931)
for(i in 1:4){
 #i <- 2
  piL <- tdList[[i]] # DIM: IT x N
  pii <- apply(piL,1,mean)  
  if(i==1){
    plot(x=median(pii),y=i, "p", pch=19, ylim=c(0.9,5), xlim=c(0.02,0.1), axes=FALSE, ylab="", xlab=expression(p[t]^TD), 
         cex.lab=1.5) # take mean of all 
    lines(x=quantile(pii, probs = c(0.025,0.975)), y=c(i,i), lwd=2)
    dd <- density(pii)
    points(x=emp[i], y=i, pch=19, col="red")
    lines(x=dd$x, y=dd$y*0.005+i, lty=1, col="grey35")
    #axis(2, at=1:5, labels=c("", "", "", "", ""), cex.axis=1.4)
    axis(1, cex.axis=1.1)
  } else {
    points(x=median(pii), y=i, pch=19)
    lines(x=quantile(pii, probs = c(0.025,0.975)), y=c(i,i),lty=1, lwd=2)
    dd <- density(pii)
    points(x=emp[i], y=i, pch=19, col="red")
    lines(x=dd$x, y=dd$y*0.005+i, lty=1, col="grey35")
  }
  cat("Time: ",i," --- Median: ", round(median(pii),3), " --- Mean: ", round(mean(pii),3), 
      " --- Quantiles: ", round(quantile(pii, probs = c(0.025,0.975)),3), " --- HPD: ", round(getHPD(pii)[[1]],3),
      " \n ----- \n")
}
text(0.03, 4.8,"B.", cex=2.7)
