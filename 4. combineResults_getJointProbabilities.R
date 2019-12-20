################################################################################
################################################################################
## BART FOR PREDICTING WAVE 5 PARTICIPATION
## Combine results
## 28.03.2019
## Sabine Zinn
################################################################################
################################################################################

rm(list=ls())
################################################################################
# I. LOAD LIBRARIES AND DATA 
################################################################################
# load workspace from model for permanent dropout
save.image("Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project\\Results\\modelIndF_est.Rdata")
# load workspace from model for temporary dropout
load("Z:\\Projects\\p000139_Methoden_Survey_Psych\\BART_Project\\Results\\partModel_est.Rdata")

################################################################################
# II. COMBINE RESULTS
################################################################################
# probability of going over to individual field: pi, i.e. staying in school context is pd=1-pi
# piList
pi.1 <- piList[[1]]
DD <- datl2[datl2$wave==1,]
iD <- which(DD$indNV==0)
pii.1 <- apply(pi.1[,iD],1,mean)
pi.2 <- piList[[2]]
DD <- datl2[datl2$wave==2,]
iD <- which(DD$indNV==0)
pii.2 <- apply(pi.2[,iD],1,mean)
pi.3 <- piList[[3]]
DD <- datl2[datl2$wave==3,]
iD <- which(DD$indNV==0)
pii.3 <- apply(pi.3[,iD],1,mean)
pi.4 <- piList[[4]]
DD <- datl2[datl2$wave==4,]
iD <- which(DD$indNV==0)
pii.4 <- apply(pi.4[,iD],1,mean)

# participation probabilites: p1, p2,p2, p4
# tdList # temp. dropout probs
pi.1 <- tdList[[1]] 
pim.1 <- apply(pi.1,1,mean)  
pi.2 <- tdList[[2]] 
pim.2 <- apply(pi.2,1,mean) 
pi.3 <- tdList[[3]] 
pim.3 <- apply(pi.3,1,mean) 
pi.4 <- tdList[[4]] 
pim.4 <- apply(pi.4,1,mean) 

# Wave 2
ps2 <- 1-median(pii.1)
p2 <- 1-median(pim.1)
pw2 <- ps2*p2
ppd.2  <- (1-pii.1)*(1-pim.1)
ddW2 <- density(ppd.2)
q.2 <- quantile(ppd.2, probs = c(0.025, 0.975))

# Wave 3
ps3 <- 1-median(pii.2)
p3 <- 1-median(pim.2)
pw3 <- ps2*p2*ps3*p3 + ps2*(1-p2)*ps3*p3
ppd.3 <- (1-pii.1)*(1-pim.1)*(1-pii.2)*(1-pim.2) + (1-pii.1)*(pim.1)*(1-pii.2)*(1-pim.2)
ddW3 <- density(ppd.3)
q.3 <- quantile(ppd.3, probs = c(0.025, 0.975))

# Wave 4
ps4 <- 1-median(pii.3)
p4 <- 1-median(pim.3)
pw4 <- ps2*p2*ps3*p3*ps4*p4 + 
        ps2*(1-p2)*ps3*p3*ps4*p4 +
          ps2*p2*ps3*(1-p3)*ps4*p4   
ppd.4 <- (1-pii.1)*(1-pim.1)*(1-pii.2)*(1-pim.2)*(1-pii.3)*(1-pim.3) + 
  (1-pii.1)*(pim.1)*(1-pii.2)*(1-pim.2)*(1-pii.3)*(1-pim.3) +
  (1-pii.1)*(1-pim.1)*(1-pii.2)*(pim.2)*(1-pii.3)*(1-pim.3)  
ddW4 <- density(ppd.4)
q.4 <- quantile(ppd.4, probs = c(0.025, 0.975))

# WAVE 5
ps5 <- 1-median(pii.4)
p5 <- 1-median(pim.4)
pw5 <- ps2*p2*ps3*p3*ps4*p4*ps5*p5 + 
  ps2*(1-p2)*ps3*p3*ps4*p4*ps5*p5 +
  ps2*p2*ps3*(1-p3)*ps4*p4*ps5*p5 +  
  ps2*p2*ps3*p3*ps4*(1-p4)*ps5*p5 +
  ps2*(1-p2)*ps3*p3*ps4*(1-p4)*ps5*p5 
ppd.5 <- (1-pii.1)*(1-pim.1)*(1-pii.2)*(1-pim.2)*(1-pii.3)*(1-pim.3)*(1-pii.4)*(1-pim.4) + 
  (1-pii.1)*(pim.1)*(1-pii.2)*(1-pim.2)*(1-pii.3)*(1-pim.3)*(1-pii.4)*(1-pim.4) +
  (1-pii.1)*(1-pim.1)*(1-pii.2)*(pim.2)*(1-pii.3)*(1-pim.3)*(1-pii.4)*(1-pim.4) +  
  (1-pii.1)*(1-pim.1)*(1-pii.2)*(1-pim.2)*(1-pii.3)*(pim.3)*(1-pii.4)*(1-pim.4) +  
  (1-pii.1)*(pim.1)*(1-pii.2)*(1-pim.2)*(1-pii.3)*(pim.3)*(1-pii.4)*(1-pim.4)  
ddW5 <- density(ppd.5)
q.5 <- quantile(ppd.5, probs = c(0.025, 0.975))

# Plot Probabilities
#W2    
plot(x=pw2,y=1, "p", pch=19, ylim=c(0.9,5), xlim=c(0.6,0.92), axes=FALSE, ylab="", 
     xlab=expression(1-p[t]^NR), cex.lab=1.5)
lines(x=as.numeric(summary(ddW2$x)[c(2,5)]), y=c(1,1), lwd=2)
lines(x=ddW2$x, y=ddW2$y*0.005+1, lty=1, col="grey35")
#axis(2, at=1:5, labels=c("Wave 2", "Wave 3", "Wave 4", "Wave 5", ""), cex.axis=1.3)
#axis(2, at=1:5, labels=c("", "", "", "", ""), cex.axis=1.3)
axis(1, cex.axis=1.1)
#W3
points(x=pw3, y=2, pch=19)
lines(x=as.numeric(summary(ddW3$x)[c(2,5)]), y=c(2,2), lwd=2)
lines(x=ddW3$x, y=ddW3$y*0.005+2, lty=1, col="grey35")
#W4
points(x=pw4, y=3, pch=19)
lines(x=as.numeric(summary(ddW4$x)[c(2,5)]), y=c(3,3), lwd=2)
lines(x=ddW4$x, y=ddW4$y*0.005+3, lty=1, col="grey35")
#W3
points(x=pw5, y=4, pch=19)
lines(x=as.numeric(summary(ddW5$x)[c(2,5)]), y=c(4,4), lwd=2)
lines(x=ddW5$x, y=ddW5$y*0.005+4, lty=1, col="grey35")
text(0.61, 4.8,"C.", cex=2.7)

