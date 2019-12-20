################################################################################
################################################################################
## BART FOR PREDICTING PROBABILITIES OF TEMPORARY DROPOUT
## Level: school context
## 18.03.2019
## Sabine Zinn
################################################################################
################################################################################


################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################
library(BART, lib.loc="Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project")
library(mice)
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

################################################################################
# II. DEFINE MODEL MATRICES AND SETTINGS
################################################################################
# Matrix of covariates as training data
x.train <- Dat[!duplicated(Dat$ID_t),c(2:84)]
colnames(x.train) <- colnames(Dat[,c(2:84)])

# Define matrix for times and events (temporary dropout), note that there are no final dropouts in the school context
ids <- unique(Dat$ID_t) # N=4559
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

# Define init structure for BART
pre <- recur.pre.bart(times=times[,-1], delta=delta[,-1], x.train= as.matrix(x.train))

# Fit model
post <- recur.bart(y.train=pre$y.train, pre$tx.train, x.test=pre$tx.test, 
                    #  ntree = 300, usequants=TRUE, nskip=200, sparse=TRUE, 
                    #  ndpost=1500, keepevery=200)
                    ntree = 300, usequants=TRUE, nskip=500, sparse=TRUE, 
                    ndpost= 5000, keepevery=500)
save.image("Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project\\Results\\partModel_est.Rdata")

################################################################################
# III. CHECK CONVERGENCE
################################################################################
# A. There is some autocorrelation: 
# (1) Close chain elements are very similar, meaning that throwing one away does not loose a lot of information (that is what the autocorrelation plot shows)
# (2) You need a lot of repetitions to get convergence, meaning that you get very large chains if you don't thin. Because of that, working with the full chain can be very slow, 
#     costs a lot of storage, or even lead to memory problems when monitoring a lot of variables.
eT <-  cbind(unique(Dat$ID_t),as.matrix(x.train))
colnames(eT)[1] <- "ID_t"
preE <- recur.pre.bart(times=times[,-1], delta=delta[,-1], x.train= eT)
K <- pre$K
allID <- preE$tx.train[preE$tx.train[,"t"]==4, "ID_t"] # choose five individuals that did not perm. dropout before study end, to compute acf
ii <- sample(allID,5)
j. <- seq(-0.5, 0.4, length.out=5)
par(mfrow = c(2,2))
for(j in 1:K) {
  for(i in c(1:5)) {
    auto.corr <- acf(post$yhat.train[ ,which(preE$tx.train[,"ID_t"]==ii[i])[j]], plot=FALSE)
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
ii <- sample(allID,100)
for(j in 1:K) {
  for(i in 1:length(allID)) {
    if(i==1)
      plot(post$yhat.train[,which(preE$tx.train[,"ID_t"]==ii[i])[j]], type='l',
           ylim=c(-3, 0), sub=paste0('t=', pre$times[j]),
           ylab=expression(Phi(g(x, T, M))), xlab='Z?ge')
    else
      lines(post$yhat.train[,which(preE$tx.train[,"ID_t"]==ii[i])[j]], type='l', col=i)
  }
  Sys.sleep(1)
}
# C. Geweke convergence plot: diagnostic for probit BART: 
LL <- 50
ii <- sample(allID,LL)
geweke <- list(1:LL)
for(i in 1:LL) {
  geweke[[i]] <- gewekediag(post$yhat.train[ ,which(preE$tx.train[,"ID_t"]==ii[i])])
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
################################################################################
# Variable importance
pdf("importExt_pM.pdf", width=25, height = 10)
par(mar=c(25, 5, 4, 2))
P <-  length(post$varprob.mean)
pchS <- ifelse(post$varprob.mean > 1/P, 19, 1)
colS <- c("black", gray(0.25), gray(0.35), gray(0.45), gray(0.55), gray(0.65))
colP <- c(rep(colS[1],3), rep(colS[2], 27), rep(colS[3],6), 
          rep(colS[4],27), rep(colS[5],7), rep(colS[6],15))
plot(post$varprob.mean, ylab='Selection Probability', ylim=c(0, 0.6),
     xlab="", pch=pchS, col=colP, axes=FALSE, cex.lab=1.2)
axis(2)
labsO <- names(post$varprob.mean)
labs <- c("t","v", "N", "sex", "age", "migration background: yes",
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
axis(1, at=1:P, labels=labs, las=2, cex.lab=0.8)
lines(c(0, length(names(post$varprob.mean))), c(1/P, 1/P))
legend("topleft", c("Wave, sojourn time, number previous dropouts",
                     "Student attributes", "Missing value dummies",
                     "Aggregrated student information on school level",
                     "School context information",
                     "Federal state dummies"), fill=colS, ncol = 1,
                      cex=1.2, bty="n")
dev.off()

pdf("importExt_red_pM.pdf", width=12, height = 10)
par(mar=c(25, 5, 4, 2))
selI <- c(1,2,3,6,39,7,40,12,45,17,50,23,56,24,57,30,63,65,67)
redI <- post$varprob.mean[selI]
P <-  length(post$varprob.mean)
pchS <- ifelse(post$varprob.mean > 1/P, 19, 1)[selI]
plot(redI, ylab='Selection Probability', ylim=c(-0.3, 0.35),
     xlab="", pch=pchS, axes=FALSE, cex.lab=1.2)
axis(2)
labsO <- names(post$varprob.mean)
labs <- c("t","v", "N", "sex", "age", "migration background: yes",
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
axis(1, at=1:length(redI), labels=labs[selI], las=2, cex.lab=0.8)
lines(c(0, length(names(redI))), c(1/P, 1/P))
dev.off()

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
RES <- cbind(colnames(XX)[ii[,2]], rownames(XX)[ii[,1]], XX[ii]) # most frequent branching after "time"
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
#names(post$varcount.mean)[81]
#post$treedraws$cutpoints$STATESH[1]
#post$treedraws$cutpoints$t[1]

################################################################################
# V. PREDICTION
################################################################################
ppO <- predict(post, pre$tx.train[,post$rm.const])   # on the observed

# Predict probabilities of temp. dropout
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
      " --- Quantiles: ", round(quantile(pii, probs = c(0.025,0.975)),3), #" --- HPD: ", round(getHPD(pii)[[1]],3),
      " \n ----- \n")
}
text(0.03, 4.8,"B.", cex=2.7)

# save all
save.image("Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project\\Results\\partModel_est.Rdata")

