##---------------------------------------------------
## Simonsen-Churchill model: Density of right tree height
## Lemma 2, equation (5) in Hobolth and Jensen (2014)
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
## Name: ARGTrns
## Purpose: Transition density for ARG
##          Density for new tree conditional on recombination
##          and left tree height
## Input: rho: recombination rate
##        Value of fixed left site tree: s
##        Values of right site coalescent smaller and
##        larger than s: t.sml and t.lrg
## Output: List:
##         Densities q.sml and q.lrg corresponding to t.sml and t.lrg
## Comments: This is eq (5) in Hobolth and Jensen (2014, TPB)
##--------------------------------------------------------
ARGTrns <- function(rho,s,t.sml,t.lrg){
  ARGRateMat <- ARGRateM(rho)
  ## s<t
  expm.s <- expm(ARGRateMat*s)
  q.lrg <- exp(-(t.lrg-s))*(expm.s[1,2]+expm.s[1,3])/
    (expm.s[1,2]+expm.s[1,3]+expm.s[1,4]+expm.s[1,6])
  ## s>t
  n <- length(t.sml)
  denom <- rep(0,n)
  for (i in 1:n){
    expm.t <- expm(ARGRateMat*t.sml[i])
    denom[i] <- expm.t[1,2]+expm.t[1,3]
  }
  q.sml <- exp(-(s-t.sml))*denom/
    (expm.s[1,2]+expm.s[1,3]+expm.s[1,4]+expm.s[1,6])
  ## Prepare output
  out <- list()
  out$q.lrg <- q.lrg
  out$q.sml <- q.sml
  return(out)
}
## Reproduce ARG part of Figure 5 in Hobolth and Jensen (2014, TPB)
s <- 0.5 ; n <- 100
t.sml <- seq(from=0.001,to=s,length=n)
t.lrg <- seq(from=s,to=4.0,length=n)
rho <- 0.02
ARGdns <- ARGTrns(rho,s,t.sml,t.lrg)
#pdf("NewTreeHeight.pdf",width=8.0,height=8.0)
plot(c(t.sml,t.lrg),c(ARGdns$q.sml,ARGdns$q.lrg),
     ylim = c(0,0.9),cex.axis=1.3,cex.lab=1.3,
     type="l",lwd=2,col="darkblue",
     xlab="Right tree height t",
     ylab="Density of right tree height")
abline(v=s,col="black")
text(0.65,0.85,expression(rho),cex=1.3,col="darkblue")
text(0.65,0.85,pos=4,"=0.02 and s=0.5",cex=1.3,col="darkblue")
#rho1 <- 0.2 ; ARGdns <- ARGTrns(rho1,s,t.sml,t.lrg)
#points(c(t.sml,t.lrg),c(ARGdns$q.sml,ARGdns$q.lrg),
#       type="l",lwd=2,col="darkblue")
#rho2 <- 0.5 ; ARGdns <- ARGTrns(rho2,s,t.sml,t.lrg)
#points(c(t.sml,t.lrg),c(ARGdns$q.sml,ARGdns$q.lrg),
#       type="l",lwd=2,col="blue")
rho3 <- 2.5 ; ARGdns <- ARGTrns(rho3,s,t.sml,t.lrg)
points(c(t.sml,t.lrg),c(ARGdns$q.sml,ARGdns$q.lrg),
       type="l",lwd=2,col="blue")
text(0.75,0.8,expression(rho),cex=1.3,col="blue")
text(0.75,0.8,pos=4,"=2.5 and s=0.5",cex=1.3,col="blue")
s <- 2 ; n <- 100
t.sml <- seq(from=0.001,to=s,length=n)
t.lrg <- seq(from=s,to=4.0,length=n)
rho <- 0.02
ARGdns <- ARGTrns(rho,s,t.sml,t.lrg)
points(c(t.sml,t.lrg),c(ARGdns$q.sml,ARGdns$q.lrg),
       type="l",lwd=2,col="darkred")
abline(v=s,col="black")
#rho1 <- 0.2 ; ARGdns <- ARGTrns(rho1,s,t.sml,t.lrg)
#points(c(t.sml,t.lrg),c(ARGdns$q.sml,ARGdns$q.lrg),
#       type="l",lwd=2,col="darkblue")
#rho2 <- 0.5 ; ARGdns <- ARGTrns(rho2,s,t.sml,t.lrg)
#points(c(t.sml,t.lrg),c(ARGdns$q.sml,ARGdns$q.lrg),
#       type="l",lwd=2,col="blue")
rho3 <- 2.5 ; ARGdns <- ARGTrns(rho3,s,t.sml,t.lrg)
points(c(t.sml,t.lrg),c(ARGdns$q.sml,ARGdns$q.lrg),
       type="l",lwd=2,col="red")
text(2.15,0.65,expression(rho),cex=1.3,col="darkred")
text(2.15,0.65,pos=4,"=0.02 and s=2",cex=1.3,col="darkred")
text(2.25,0.6,expression(rho),cex=1.3,col="red")
text(2.25,0.6,pos=4,"=2.5 and s=2",cex=1.3,col="red")
#dev.off()