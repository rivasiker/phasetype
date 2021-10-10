# Name: logFrwdLikFct
# Written: 8 March 2006
# Author: Asger Hobolth
# Purpose:
# Calculates the likelihood function using forward recursion
#--------------------------------------------------------
# Input:
# InitProb: Initial probabilities:
#           Probabilities of going from the start state
#           The start state does not emit any letters.
# TransProb: Transition probabilities;
#            Probabilities between hidden states
#            nHS times nHS matrix
# EndProb: Ending probabilities
#          Probabilities of going to the end state
#          The end state does not emit any letters.
# EmisProb: Emission probabilities
# ObsSeq: Observed sequence
#
# Output:
# logForwardLikVal: log likelihood value
#---------------------------------------------------------
"logFrwdLikFct" <- function(InitProb,TransProb,EndProb,EmisProb,ObsSeq){
  len <- length(ObsSeq)
  nHS <- nrow(TransProb)
  #----------------------------------
  # Forward algorithm with log values
  #----------------------------------
  # Define logForwardLik
  logForwardLik <- matrix(0,nrow=len,ncol=nHS)
  # Start condition
  logForwardLik[1,] <- log(InitProb*EmisProb[,ObsSeq[1]])
  # Determine logForwardLik by recursion
  for (k in 2:len){
    a <- max(logForwardLik[k-1,])
    for (j in 1:nHS){
      logForwardLik[k,j] <- log(sum(TransProb[,j]*
                                      rep(EmisProb[j,ObsSeq[k]],nHS)*
                                      exp(logForwardLik[k-1,]-a)))+a
    }
  }
  a <- max(logForwardLik[len,])
  logForwardLikVal <- log(sum(EndProb*exp(logForwardLik[len,]-a)))+a
  cat("log likelihood from forward algorithm:",logForwardLikVal,"\n")
  return(logForwardLikVal)
}