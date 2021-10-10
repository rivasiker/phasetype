source("HMMlogLik.R")
##-----------------------------------------------------------------
## Main 
##-----------------------------------------------------------------
#setwd("C:/Users/asger/Projects/ParticleFilter/Programs/")
## Need expm package: install.packages("expm")
#library("expm")
## Define time intervals, recombination rate and transition matrix
nInt <- 30 ; tm <- -log(1-1:(nInt)/nInt)
## Need transition matrix
source("TransMat.R")
rho <- 0.1 ; TransM <- FastTransMat(tm,rho)
#rho <- 2.5 ; TransM <- FastTransMat(tm,rho)
##---------------------------------------------
## IKER: WE DISCUSSED UP TO THIS POINT
##---------------------------------------------
## Sequence length and tht
sq.len <- 20000 ; tht <- 0.1
source("Sim.R")
sim <- SimFct(TransM,tm,tht,sq.len)
##------------------------------------------------------------------
## Plot tree height
plot(tm[sim$hght.sq],ylim=c(0,4),lwd=2,type="s",
     xlab="Genomic position",ylab="Discretized time")
## Add mutations to the plot
het.indx <- which(sim$dipl.sq==1)
n.hets <- length(het.indx)
points(het.indx,rep(0,n.hets),pch=19,cex=0.2,col="red")
##------------------------------------------------------------------
## Calculate Ripley's empirical K function 
## (actually, probability of mutation at a certain distance)
RipleysEmpirFct <- function(mx,het.indx){
  #mx <- 50
  tbl <- rep(0,mx)
  tbl.add <- 1 ; i <- 0
  while (sum(tbl.add)!=0){
    i <- i+1
    tbl.add <- tabulate(diff(het.indx,i),nbins=mx)
    tbl <- tbl+tbl.add
    #cat("Adds in tbl:",sum(tbl.add),"\n")
    #cat("i:",i,"\n")
    if (i>40) break
  }
  ## Observation vector
  emp.v <- tbl
  return(tbl)
}
obs.v <- RipleysEmpirFct(50,het.indx)
source("RipleysK.R")
source("StDst.R")
source("RipleysKfitRho.R")
## Estimate rho using Ripleys function
Ripley.rho.est <- FitRhoFct(obs.v,length(het.indx))
## Plot empirical and expected Ripleys function 
TransM <- TransMat(tm,Ripley.rho.est)
mx <- 50
expc.prb <- RipleysMutafct(TransM,tm,tht)$muta.prb[1:mx]
plot(1:mx,expc.prb,type="l",col="blue",
     main="Probability of mutation",xlab="Distance",ylab="Probability")
points(1:mx,obs.v/length(het.indx),pch=19,cex=0.5,col="red")
##------------
## Runs of homozygosity
## Empirical distribution
mx.rh <- 10
obs.rh <- tabulate(diff(het.indx))[1:mx.rh]/(length(het.indx)-1)
source("RunsOfHomoz.R")
source("RunsOfHomozFitRho.R")
## Estimate rho using runs of homozygosity
RunsOfHomoz.rho.est <- RunsOfHomozygosityFitRhoFct(obs.rh)
## Plot empirical and expected runs of homozygosity 
plot(1:mx.rh,obs.rh,ylim=c(0,0.2))
TransM <- FastTransMat(tm,Ripley.rho.est)
expc.rh <- RunsOfHomozPrb(TransM,tm,tht)$wait.dst[1:mx.rh]
points(1:mx.rh,expc.rh,type="l",col="red")
##--------------------
## Simulation study
##--------------------
#n.sim <- 50
n.sim <- 5
res.sim <- matrix(0,nrow=n.sim,ncol=3)
for (j in 1:n.sim){
  rho <- 0.1 ; TransM <- FastTransMat(tm,rho)
  sim <- SimFct(TransM,tm,tht,sq.len)
  het.indx <- which(sim$dipl.sq==1)
  obs.v <- RipleysEmpirFct(50,het.indx)
  Ripley.rho.est <- FitRhoFct(obs.v,length(het.indx))
  mx.rh <- 10
  obs.rh <- tabulate(diff(het.indx))[1:mx.rh]/(length(het.indx)-1)
  RunsOfHomoz.rho.est <- RunsOfHomozygosityFitRhoFct(obs.rh)
  ##-------------------------------------------------------
  ## HMM Estimation
  rho.v <- seq(from=0.07,to=0.13,by=0.01)
  logL.v <- rep(0,len=length(rho.v))
  ## Loop over possible values of rho
  for (i in 1:length(rho.v)){
    rho <- rho.v[i]
    TransM <- FastTransMat(tm,rho)
    #EmisPrb <- matrix(c(exp(-tm*tht),1-exp(-tm*tht)),ncol=2,byrow=FALSE)
    lw.tm <- c(0,tm[1:(nInt-1)]) ; df.tm <- diff(c(0,tm))  
    muta.tm <- 1-exp(-tht*lw.tm)*(1-exp(-(1+tht)*df.tm))/(1-exp(-df.tm))/(1+tht)
    EmisPrb <- matrix(c(1-muta.tm,muta.tm),ncol=2,byrow=FALSE)
    ## log Likelihood
    logL.v[i] <- 
      logFrwdLikFct((1:nInt)/nInt,TransM,(1:nInt)/nInt,EmisPrb,(1+sim$dipl.sq))
  }
  ## Estimate HMM rho 
  spl <- smooth.spline(rho.v,logL.v)
  pred.spl <- predict(spl,x <-seq(0.07,0.13,len=200))
  plot(rho.v,logL.v)
  points(pred.spl$x,pred.spl$y,type="l",col="red")
  HMM.rho.est <- pred.spl$x[which.max(pred.spl$y)]
  cat("Sim:",j,"Estimates:",Ripley.rho.est,RunsOfHomoz.rho.est,HMM.rho.est,"\n")
  res.sim[j,] <- c(Ripley.rho.est,RunsOfHomoz.rho.est,HMM.rho.est)
}
plot(1:n.sim,sort(res.sim[,1]),ylim=c(0.005,0.195),pch=19,cex=.5,
     main="Simulation study",xlab="Simulation Number",ylab="Estimate of Recombination Rate")
points(1:n.sim,sort(res.sim[,1]),type="l")
abline(h=0.1)
points(1:n.sim,pmin(sort(res.sim[,2]),0.195),col="red",pch=19,cex=.5)
points(1:n.sim,pmin(sort(res.sim[,2]),0.195),col="red",type="l")
points(1:n.sim,sort(res.sim[,3]),col="blue",pch=19,cex=.5)
points(1:n.sim,sort(res.sim[,3]),col="blue",type="l")
legend("topleft",c("Hidden Markov Model","Pair Correlation","Nearest Neighbour"),
       col=c("blue","black","red"),lty=1,bty="n")
#write(res.sim,file="ResSim.dat",ncolumns=3)
