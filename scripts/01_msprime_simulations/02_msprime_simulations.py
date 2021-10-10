#!/usr/bin/env python
# coding: utf-8

# In[1]:


import msprime
import numpy as np
import sys
from IPython.display import SVG
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
# %reload_ext rpy2.ipython


# In[6]:


ide = sys.argv[1]
rho = sys.argv[2]
tht = sys.argv[3]


# In[7]:


get_ipython().run_cell_magic('R', '-i rho -i tht -i ide', '\nlibrary(tidyverse)\nlibrary(expm)\nlibrary(parallel)\nlibrary(optimParallel)')


# In[8]:


ide = sys.argv[1]
rho = float(sys.argv[2])
tht = float(sys.argv[3])


# In[9]:


get_ipython().run_cell_magic('R', '', "\nFastTransMat <- function(tm,rho){\n  # This produces a rate matrix which has an additional state,\n  # number 4, which represents states 4, 5, 6 and 7 of the \n  # 'slow' matrix. Thus, state 5 corresponds to the absorbing\n  # state, i.e. state 8 in the 'slow' matrix. \n  rate_mat <- matrix(c(-(1+rho),rho,0,0,1,\n                              1,-(3+rho/2),rho/2,2,0,             \n                              0,4,-6,2,0,\n                              0,0,0,-1,1,\n                              0,0,0,0,0),\n                            nrow=5,ncol=5,byrow=TRUE)\n  nInt <- length(tm) ## Number of intervals\n  tm0 <- c(0,tm) ## tm0 is tm with time 0 added\n  ##-------------------------------------\n  JointMat <- matrix(0,nrow=nInt,ncol=nInt) ## Joint prb matrix\n  for (j in 1:(nInt-1)){  ## Left state\n    for (k in j:nInt){  ## Right state\n      if (j<k){\n        JointMat[j,k] <-\n          expm(tm0[j]*rate_mat)[1,1:3]%*%\n          expm((tm[j]-tm0[j])*rate_mat)[1:3,4]*\n          exp(-(tm0[k]-tm[j]))*\n          # which is the same as \n          # expm((tm0[k]-tm[j])*rate_mat)[4,4]*0.5*\n            (1-exp(-(tm[k]-tm0[k])))*0.5\n            # which is the same as \n            # expm((tm[k]-tm0[k])*rate_mat)[4,5]\n            # but because expm can't handle infinity, \n            # we need to do a workaround \n        # Symmetrize\n        JointMat[k,j] <- JointMat[j,k]\n      }\n      if (j==k){\n        JointMat[j,k] <-\n          expm(tm0[j]*rate_mat)[1,1:3]%*%\n          expm((tm[j]-tm0[j])*rate_mat)[1:3,5]\n      }\n    }\n  }\n  ## Final entry\n  JointMat[nInt,nInt] <- sum(expm(tm0[nInt]*rate_mat)[1,1:3])\n  ## Again: expm can't handle infinity, and state 8 is absorbing\n  ## Transition matrix \n  TransMat <- JointMat/rowSums(JointMat)\n  return(TransMat)\n}")


# In[10]:


get_ipython().run_cell_magic('R', '', '\ncalc_p_mut <- function(tm, tht) {\n    tm_1 <- c(0, tm[1:(length(tm)-1)])\n    prob_mutation <- rep(NA, length(tm))\n    for (i in 1:length(tm)){\n        prob_mutation[i] <- 1-exp(-tht*tm_1[i])*(1-exp(-(1+tht)*(tm[i]-tm_1[i])))/(1-exp(-(tm[i]-tm_1[i])))/(1+tht)\n    }\n    prob_mutation\n}')


# In[11]:


get_ipython().run_cell_magic('R', '', '\nlogFrwdLikFct <- function(InitProb,TransProb,EndProb,EmisProb,ObsSeq){\n  # Number of observations\n  len <- length(ObsSeq)\n  # Number of hidden states\n  nHS <- nrow(TransProb)\n  # Define logForwardLik\n  logForwardLik <- matrix(0,nrow=len,ncol=nHS)\n  # Start condition\n  logForwardLik[1,] <- log(InitProb*EmisProb[,ObsSeq[1]])\n  # Determine logForwardLik by recursion \n  for(k in 2:len){\n    a <- max(logForwardLik[k-1,])\n    # Calculate the loglik of current iteration\n    logForwardLik[k,] <- \n      log(colSums(\n          (exp(logForwardLik[k-1,]-a)%*%TransProb)*EmisProb[,ObsSeq[k]]\n      ))+a\n  }  \n  a <- max(logForwardLik[len,])\n  logForwardLikVal <- log(sum(EndProb*exp(logForwardLik[len,]-a)))+a\n  return(logForwardLikVal)\n}')


# In[12]:


ts = msprime.sim_ancestry(
    samples=2,
    recombination_rate=rho/2,
    sequence_length=100_000,
    population_size = 1,
    ploidy = 1
)

mutated_ts = msprime.sim_mutations(
    ts,
    rate = tht/2,
    model = 'infinite_alleles'
)


# In[13]:


tree_heights = ts.tables.nodes[2:].asdict()['time'].tolist()


# In[14]:


mut_mat = np.zeros((mutated_ts.num_nodes-2, 2))
first = True
# For each tree
for tree in mutated_ts.trees():
    # Save the number of positions in the tree
    mut_mat[tree.root-2, 0] += tree.span
    # save the number of mutations
    mut_mat[tree.root-2, 1] += tree.num_mutations


# In[30]:


mut_pos = []
for variant in mutated_ts.variants():
    mut_pos.append(int(variant.position))
mut_pos_bin = np.zeros(int(mutated_ts.sequence_length), dtype=np.int32)
for i in range(int(mutated_ts.sequence_length)):
    if i in mut_pos:
        mut_pos_bin[i] = 1
    else:
        mut_pos_bin[i] = 0


# In[39]:


get_ipython().run_cell_magic('R', '-i mut_pos_bin', "\ncalc_tm <- function(nInt) -log(1-1:(nInt)/nInt)\nnInt <- 20\ntm <- calc_tm(nInt)\n\n\n\nfun_optim <- function(x) {\n    require(expm)\n    rho <- x[1]\n    tht <- x[2]\n    initial <- rep(1/nInt, nInt)\n    transition <- FastTransMat(tm, rho)\n    emission <- matrix(c(1-calc_p_mut(tm, tht),calc_p_mut(tm, tht)),ncol=2,byrow=FALSE)\n    logFrwdLikFct(initial, transition, initial, emission, mut_pos_bin+1)\n}\n\ncl <- makeCluster(20)     # set the number of processor cores\nsetDefaultCluster(cl=cl)  # set 'cl' as default cluster\nclusterExport(cl=cl, c('nInt', 'mut_pos_bin', 'tm', 'FastTransMat', 'calc_p_mut', 'logFrwdLikFct'))\noptim_both <- optimParallel(c(0.1, 0.01), \n                            fun_optim, \n                            lower = c(0.000001, 0.000001),\n                            upper = c(10, 10),\n                            parallel = list(loginfo = TRUE),\n                            control = list(ndeps=1e-4, fnscale = -1, pgtol = 0),\n                            )")


# In[40]:


get_ipython().run_cell_magic('R', '', '\noptim_both')


# In[19]:


get_ipython().run_cell_magic('R', '', "\nwrite_csv(as_tibble(optim_both$loginfo), paste0('../../steps/01_msprime_simulations/sim', ide, '_rho', rho, '_tht', tht, '.csv'))")


# In[1]:


get_ipython().system('jupyter nbconvert --to script 02_msprime_simulations.ipynb')


# In[ ]:




