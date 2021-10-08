#setwd("C:/Users/asger/Projects/ParticleFilter/Programs/")
##---------------------------------------------------
## Simonsen-Churchill model: Probability of a change
## Lemma 2, equation (4) in Hobolth and Jensen (2014)
##---------------------------------------------------
## Need expm package: install.packages("expm")
library("expm")
##--------------------------------------------------------
##----------- Simonsen-Churchill rate matrix -------------
##--------------------------------------------------------
## Name: ARGRateM
## Input: rho: recombination rate
## Output: 8x8 Simonsen-Churchill rate matrix
## Comments: This is the rate matrix in Fig 1a in
##           Hobolth and Jensen (2014, TPB)
##--------------------------------------------------------
ARGRateM <- function(rho){
  RateM <-
    matrix(c(0,rho,0,0,0,0,0,1,
             1,0,rho/2,1,1,0,0,0,             
             0,4,0,0,0,1,1,0,
             0,0,0,0,0,rho/2,0,1,
             0,0,0,0,0,0,rho/2,1,
             0,0,0,2,0,0,0,1,
             0,0,0,0,2,0,0,1,
             0,0,0,0,0,0,0,1,
             0,0,0,0,0,0,0,0),nrow=8,ncol=8,byrow=TRUE)
  ## Get diagonals right
  for (rw in 1:8){
    RateM[rw,rw] <- -sum(RateM[rw,])
  }
  return(RateM)
}
## Example
rho <- 0.1 ; ARGRateMat <- ARGRateM(rho)
##--------------------------------------------------------
## Name: ARGHoldingTime
## Purpose: Holding time for ARG (prb for staying in state)
## Input: Parameter: rho
##        Values of left locus coalescent: s.vec
## Output: Probability vector p.vec corresponding to s.vec
## Comments: This is eq (4) in Hobolth and Jensen (2014, TPB)
##--------------------------------------------------------
ARGHoldingTime <- function(rho,s.vec){
  ARGRateMat <- ARGRateM(rho)
  n <- length(s.vec)
  p.vec <- rep(0,n)
  for (i in 1:n){
    p.vec[i] <- expm((ARGRateMat*s.vec[i]))[1,1]/exp(-s.vec[i])
  }
  return(p.vec)
}
## Reproduce ARG part of Figure 4 in Hobolth and Jensen (2014, TPB)
n <- 100 ; s.vec <- seq(from=0.001,to=4.0,length=n) 
rho <- 0.02 ; p.vec <- ARGHoldingTime(rho,s.vec)
rho1 <- 0.2 ; p.vec1 <- ARGHoldingTime(rho1,s.vec)
rho2 <- 0.5 ; p.vec2 <- ARGHoldingTime(rho2,s.vec)
rho3 <- 2.5 ; p.vec3 <- ARGHoldingTime(rho3,s.vec)
#pdf("ChangePrb.pdf",width=8.0,height=8.0)
plot(s.vec,p.vec,
     xlab="Left tree height s",ylab="Probability of same tree height",ylim=c(0,1),
     type="l",lwd=2,col="black",
     cex.axis=1.3,cex.lab=1.3)
points(s.vec,p.vec1,col="darkblue",type="l",lwd=2)
points(s.vec,p.vec2,col="blue",type="l",lwd=2)
points(s.vec,p.vec3,col="red",type="l",lwd=2)
text(0,0.4,expression(rho),cex=1.3,col="red") 
text(0,0.4,paste("=",round(rho3,digits=2)),pos=4,cex=1.3,col="red")
text(1,0.6,expression(rho),cex=1.3,col="blue")
text(1,0.6,paste("=",round(rho2,digits=2)),pos=4,cex=1.3,col="blue")
text(3,0.65,expression(rho),cex=1.3,col="darkblue")
text(3,0.65,paste("=",round(rho1,digits=2)),pos=4,cex=1.3,col="darkblue")
text(3.6,0.92,expression(rho),cex=1.3,col="black")
text(3.6,0.92,paste("=",round(rho,digits=2)),pos=4,cex=1.3,col="black")
#dev.off()