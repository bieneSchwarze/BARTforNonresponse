################################################################################
################################################################################
## BART FOR PREDICTING PARTICIPATION PROBABILITIES
## Level: school context
## 18.03.2019
## Sabine Zinn
################################################################################
################################################################################

################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################
rm(list=ls())
set.seed(654)
library(BART)
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
x.train <- D[!duplicated(D$ID_t), c(2:84)]
colnames(x.train) <- colnames(D[,c(2:84)]) 

ids <- unique(D$ID_t) # N=4559
N <- length(ids)
L <- 4
times <- rep(NA, N)
delta <- rep(NA, N)
for(i in 1:length(ids)){
  id <- ids[i]
  mm <- D[D$ID_t %in% id,]
  if( 1 %in% mm$indNV){
    times[i] <- mm$time[which(mm$indNV==1)]
    delta[i] <- 1
  } else{
    #cat("It: ",i,"\n")
    times[i] <- mm$time[nrow(mm)]
    delta[i] <- 0
  }
}
table(times,delta)

post <- surv.bart(x.train, times=times, delta=delta,
                  ntree = 300, usequants=TRUE, nskip=500, sparse=TRUE, 
                  ndpost=2000, keepevery=500)
save.image("Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project\\Results\\modelIndF_est.Rdata")

################################################################################
# III. CHECK CONVERGENCE
################################################################################
# A. There is some autocorrelation: 
# (1) Close chain elements are very similar, meaning that throwing one away does not loose a lot of information (that is what the autocorrelation plot shows)
# (2) You need a lot of repetitions to get convergence, meaning that you get very large chains if you don't thin. Because of that, working with the full chain can be very slow, 
#     costs a lot of storage, or even lead to memory problems when monitoring a lot of variables.
pre <- surv.pre.bart(times=times, delta=delta, x.train=as.matrix(x.train), x.test=as.matrix(x.train))
pred <- predict(post, pre$tx.test[,-1])
k <- floor(seq(1, N, length.out=5))
j. <- seq(-0.5, 0.4, length.out=5)
K <- pre$K
par(mfrow = c(2,2))
for(j in 1:K) {
  for(i in 1:5) {
    h <- (k[i]-1)*K+j
    auto.corr <- acf(pred$yhat.test[ , h], plot=FALSE)
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
LL <- 100
k <- floor(seq(1, N, length.out=LL))
for(j in 1:K) {
  for(i in 1:LL) {
    h <- (k[i]-1)*K+j
    if(i==1)
      plot(pred$yhat.test[ , h], type='l',
           ylim=c(-4, 4), sub=paste0('t=', pre$times[j]),
           ylab=expression(Phi(g(x,T,M))), xlab='m')
    else
      lines(pred$yhat.test[ , h], type='l', col=i)
  }
  Sys.sleep(1)
}
# C. Geweke convergence plot: diagnostic for survival model with BART
LL <- 50
k <- floor(seq(1, N, length.out=LL))
K <- pre$K
geweke <- list(1:LL)
for(i in 1:LL) {
  h <- (k[i]-1)*K+1:K
  geweke[[i]] <- gewekediag(pred$yhat.test[ ,h])
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
    text(c(0.1, 0.1), c(-1.96, 1.96), pos=2, cex=0.6, labels='0.95')
    text(c(0.1, 0.1), c(-2.576, 2.576), pos=2, cex=0.6, labels='0.99')
    text(c(0.1, 0.1), c(-3.291, 3.291), pos=2, cex=0.6, labels='0.999')
    text(c(0.1, 0.1), c(-3.891, 3.891), pos=2, cex=0.6, labels='0.9999')
    text(c(0.1, 0.1), c(-4.417, 4.417), pos=2, cex=0.6, labels='0.99999')
  }
  else points(pre$times, geweke[[i]]$z)#type='l')
}
################################################################################
# IV. VARIABLE SELECTION & TREES 
################################################################################
# Variable importance
par(mar=c(7, 5, 4, 2))
P <-  length(post$varprob.mean)
pchS <- ifelse(post$varprob.mean > 1/P, 19, 1)
colS <- c("black", "red", "green", "blue", "orange", "brown")
colP <- c(rep(colS[1],1), rep(colS[2], 27), rep(colS[3],6), 
          rep(colS[4],27), rep(colS[5],7), rep(colS[6],15))
plot(post$varprob.mean, ylab='Selection Probability', ylim=c(0, 0.35),
     xlab="", pch=pchS, col=colP, axes=FALSE)
axis(2)
axis(1, at=1:P, labels=names(post$varprob.mean), las=2)
lines(c(0, length(names(post$varprob.mean))), c(1/P, 1/P))
post$varprob.mean>1/P
legend("topright", c("Wave                                   ",
                    "Student attributes", "Missing value dummies",
                    "Aggregrated student information on school level",
                    "School context information",
                    "Federal state dummies"), fill=colS, ncol=1, 
                    cex=0.9, bty="n")

# Trees
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
RES <- cbind(colnames(XX)[ii[,2]], rownames(XX)[ii[,1]], XX[ii]) # most frequent branching after covariate "time"
RES[order(RES[,3], RES[,1], decreasing = TRUE),]
colnames(RES) <- c("firstNode", "followingNode", "relFreqOfallBranches")

table(is.na(trees[-1,2])) # number of trees in total (num interations x times of trees per iteration)
tab <- table(trees[is.na(trees[,2]), 1]) # distribution of total number of nodes in all trees 
round(tab / sum(tab),3)

tab <- table(trees[which(is.na(trees[,2]))+1,2]) # freq. distr. of first node variable
names(tab) <- names(post$varcount.mean)
tab <- sort(tab, decreasing=TRUE)
round(tab / sum(tab),3)

#head(trees)
#names(post$varcount.mean)[61]
#post$treedraws$cutpoints$t[1]
#post$treedraws$cutpoints$t[1]

################################################################################
# V. PREDICTION: Probability to stay in school context, i.e. no transition in ind. field
################################################################################
pred <- predict(post, newdata=pre$tx.test[,post$rm.const])
M <- length(times)
K <- pre$K
sv <- matrix(nrow=M, ncol=K)
for(j in 1:K) {
  h <- seq(from=j, to=ncol(pred$surv.test), by=K)
  sv[ , j] <- apply(pred$surv.test[ , h], 2, mean)
}

# Plot survival curve
sv.mu  <- apply(sv, 2, mean)
diff(c(1,sv.mu))
sv.025 <- apply(sv, 2, quantile, probs=0.025)
sv.975 <- apply(sv, 2, quantile, probs=0.975)
plot(c(0, pre$times), c(1, sv.mu), type='s', ylim=0:1, ylab='WK in Schulkontext zu bleiben', xlab="", axes=FALSE)
axis(1, at=0:4, labels=c("Wave 1", "Wave 2", "Wave 3", "Wave 4", "Wave 5"))
axis(2)
lines(c(0, pre$times), c(1, sv.025), type='s', lty=2)
lines(c(0, pre$times), c(1, sv.975), type='s', lty=2)

# Derive posterior draws of prob. of temp. dropout, for each wave separately
pi.1 <- pred$prob.test[,seq(from=1, to=ncol(pred$prob.test), by=pred$K)]
hist(pi.1)
median(pi.1); mean(pi.1)
pi.2 <- pred$prob.test[,seq(from=2, to=ncol(pred$prob.test), by=pred$K)] 
hist(pi.2)
median(pi.2); mean(pi.2)
pi.3 <- pred$prob.test[,seq(from=3, to=ncol(pred$prob.test), by=pred$K)] 
hist(pi.3)
median(pi.3); mean(pi.3)
pi.4 <- pred$prob.test[,seq(from=4, to=ncol(pred$prob.test), by=pred$K)] 
hist(pi.4)
median(pi.4); mean(pi.4)
piList <- list(p1=pi.1, p2=pi.2, p3=pi.3, p4=pi.4)

# Function to derive the highest posterior density region 
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

# Have a look at posterior draws for prob. for perm. dropout for ten indiv. at wave 2
par(mfrow=c(4,3))
for(ind in 1:10){ # for all individuals (here n=10) -> one distribution for the estimated prob (i.e. one distribution per ind.)
  dd <- density(piL[,ind])
  plot(x=dd$x, y=dd$y, lty=1, type="l", col="grey35")
  points(median(piL[,ind]),0)
}

# Derive prob. of permanent dropout, averaged over all individuals
par(mfrow=c(1,3))
emp <- c(0.064, 0.108, 0.091, 0.099)
for(i in 1:4){
  #i <- 2
  piL <- piList[[i]]
  DD <- datl2[datl2$wave==i,]
  iD <- which(DD$indNV==0)
  pii <- apply(piL[,iD],1,mean)
  if(i==1){
    plot(x=median(pii),y=i, "p", pch=19, ylim=c(0.9,5), xlim=c(0,0.15), 
         axes=FALSE, ylab="", xlab=expression(p[t]^PD),cex.lab=1.5)
    # take mean of all 
    lines(x=quantile(pii, probs = c(0.025,0.975)), y=c(i,i), lwd=2)
    dd <- density(pii)
    points(x=emp[i], y=i, pch=19, col="red")
    lines(x=dd$x, y=dd$y*0.005+i, lty=1, col="grey35")
    axis(2, at=1:5, labels=c("Wave 2", "Wave 3", "Wave 4", "Wave 5", ""), cex.axis=1.3)
    axis(1, cex.axis=1.1)
  } else {
    points(x=median(pii), y=i, pch=19)
    lines(x=quantile(pii, probs = c(0.025,0.975)), y=c(i,i),lty=1, lwd=2)
    dd <- density(pii)
    points(x=emp[i], y=i, pch=19, col="red")
    lines(x=dd$x, y=dd$y*0.005+i, lty=1, col="grey35")
  }
  cat("Time: ",i," --- Median: ", round(median(pii),3), " --- Mean: ", round(mean(pii),3), 
      " --- Quantiles: ", round(quantile(pii, probs = c(0.025,0.975)),3), #" --- HPD: ", round(getHPD(pii)[[1]],3),
      " \n ----- \n")
}
text(0.03, 4.8,"A.", cex=2.7)

# save all
save.image("Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project\\Results\\modelIndF_est.Rdata")

