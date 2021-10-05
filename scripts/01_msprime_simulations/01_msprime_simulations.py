#!/usr/bin/env python
# coding: utf-8

# In[1]:


import msprime
import numpy as np
from IPython.display import SVG
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')
# %reload_ext rpy2.ipython


# # A state space model for mutations along two genomes

# ## State space model: Simonsen-Churchill framework

# This is the full transition rate matrix:

# In[2]:


get_ipython().run_cell_magic('R', '', '\nARGRateM_slow <- function(rho){\n  RateM <-\n    matrix(c(0,rho,0,0,0,0,0,1,\n             1,0,rho/2,1,1,0,0,0,             \n             0,4,0,0,0,1,1,0,\n             0,0,0,0,0,rho/2,0,1,\n             0,0,0,0,0,0,rho/2,1,\n             0,0,0,2,0,0,0,1,\n             0,0,0,0,2,0,0,1,\n             0,0,0,0,0,0,0,0),nrow=8,ncol=8,byrow=TRUE)\n  ## Get diagonals right (must sum to 0)\n  for (rw in 1:8){\n    RateM[rw,rw] <- -sum(RateM[rw,])\n  }\n  return(RateM)\n}\n\nARGRateMat_slow <- ARGRateM_slow(0.1*2)\n\nARGRateMat_slow')


# This is the fast version of the function above, where some of the states are removed from the chain if we are only interested in the heights:

# In[3]:


get_ipython().run_cell_magic('R', '', '\nARGRateM <- function(rho){\n  RateM <-\n    matrix(\n        c(-(1+rho),        rho,     0, 1,\n                 1, -(3+rho/2), rho/2, 2,\n                 0,          4,    -6, 2, \n                 0,          0,     0, 0),\n        nrow=4,ncol=4,byrow=TRUE)\n}\n\nARGRateMat <- ARGRateM(0.1*2)\n\nARGRateMat')


# Equation (1) and figure 2A can be coded this way:

# In[4]:


get_ipython().run_cell_magic('R', '', '\n\nlibrary("tidyverse")\nlibrary(\'expm\')\nlibrary(\'reshape2\')\nlibrary(\'patchwork\')')


# In[5]:


get_ipython().run_cell_magic('R', '-w 800 -h 500', "\n\n\nn <- 1000\ns <- seq(from=0.001, to=4.0,length=n)\n\n# Table with the values\ntab_dat <- expand.grid(\n        s_tib = s, \n        rho_tib = c(0.02, 0.2, 0.5, 2.5)) %>%\n    # For each row\n    rowwise() %>%\n    # Get the probablity by applying eq. (1)\n    mutate(\n        den = (exp(s_tib)*expm(ARGRateM(rho_tib)*s_tib)[1,1])) %>%\n    ungroup() %>%\n    mutate(\n        Recombination = as.character(rho_tib)\n    )\n    \n    \ntab_dat %>% \n    ggplot() + \n    geom_line(aes(s_tib, den, color = Recombination)) +\n    theme_minimal() +\n    ylab('Probability of the same tree height') +\n    xlab('Left tree height (s)')\n    ")


# The following code is for reproducing figure 2B from equation (2):

# In[6]:


get_ipython().run_cell_magic('R', '-w 800 -h 500', "\n# This function is the variable part\n# within equation (2), this is, the \n# numerator\nARGTrns <- function(rho, s, t){\n    ARGRateMat <- ARGRateM(rho)\n    expm.s <- expm(ARGRateMat*s)\n    exp(-(t-s))*(expm.s[1,2]+expm.s[1,3])\n    \n}\n\n# Define the data frame with array of t\n# and different values of rho and base s\nn <- 1000\nt <- seq(from=0.001, to=4.0,length=n)  \ntab_dat <- expand.grid(\n        t_tib = t, \n        rho_tib = c(0.02, 2.5), \n        s_tib = c(0.5, 2)) %>% \n    # For each row\n    rowwise() %>%\n    mutate(\n        # Calculate numerator\n        den = \n            ifelse(\n                t_tib > s_tib, \n                ARGTrns(rho_tib, s_tib, t_tib),\n                ARGTrns(rho_tib, t_tib, s_tib)),\n        # Divide by denominator\n        den = \n            den/(exp(-s_tib)-expm(ARGRateM(rho_tib)*s_tib)[1,1])) %>%\n    ungroup() %>%\n    mutate(\n        Recombination = as.character(rho_tib),\n        s = as.character(s_tib)\n    )\n    \ntab_dat %>% \n    ggplot() + \n    geom_line(aes(t_tib, den, color = Recombination, linetype = s)) +\n    theme_minimal() +\n    ylab('Density') +\n    xlab('Right tree height (t)') +\n    ggtitle('Density of t conditional on being different from s')\n")


# 

# ## Time discretization: setting up the finite state HMM

# In[7]:


get_ipython().run_cell_magic('R', '', "\nFastTransMat <- function(tm,rho){\n  # This produces a rate matrix which has an additional state,\n  # number 4, which represents states 4, 5, 6 and 7 of the \n  # 'slow' matrix. Thus, state 5 corresponds to the absorbing\n  # state, i.e. state 8 in the 'slow' matrix. \n  rate_mat <- matrix(c(-(1+rho),rho,0,0,1,\n                              1,-(3+rho/2),rho/2,2,0,             \n                              0,4,-6,2,0,\n                              0,0,0,-1,1,\n                              0,0,0,0,0),\n                            nrow=5,ncol=5,byrow=TRUE)\n  nInt <- length(tm) ## Number of intervals\n  tm0 <- c(0,tm) ## tm0 is tm with time 0 added\n  ##-------------------------------------\n  JointMat <- matrix(0,nrow=nInt,ncol=nInt) ## Joint prb matrix\n  for (j in 1:(nInt-1)){  ## Left state\n    for (k in j:nInt){  ## Right state\n      if (j<k){\n        JointMat[j,k] <-\n          expm(tm0[j]*rate_mat)[1,1:3]%*%\n          expm((tm[j]-tm0[j])*rate_mat)[1:3,4]*\n          exp(-(tm0[k]-tm[j]))*\n          # which is the same as \n          # expm((tm0[k]-tm[j])*rate_mat)[4,4]*0.5*\n            (1-exp(-(tm[k]-tm0[k])))*0.5\n            # which is the same as \n            # expm((tm[k]-tm0[k])*rate_mat)[4,5]\n            # but because expm can't handle infinity, \n            # we need to do a workaround \n        # Symmetrize\n        JointMat[k,j] <- JointMat[j,k]\n      }\n      if (j==k){\n        JointMat[j,k] <-\n          expm(tm0[j]*rate_mat)[1,1:3]%*%\n          expm((tm[j]-tm0[j])*rate_mat)[1:3,5]\n      }\n    }\n  }\n  ## Final entry\n  JointMat[nInt,nInt] <- sum(expm(tm0[nInt]*rate_mat)[1,1:3])\n  ## Again: expm can't handle infinity, and state 8 is absorbing\n  ## Transition matrix \n  TransMat <- JointMat/rowSums(JointMat)\n  return(TransMat)\n}")


# In[8]:


get_ipython().run_cell_magic('R', '', '\n\ncalc_tm <- function(nInt) -log(1-1:(nInt)/nInt)\n\nnInt <- 20\ntm <- calc_tm(nInt)\nrho <- 0.1*2\nTransM <- FastTransMat(tm, rho)\n\nprint(tm)\nprint(diag(TransM))\nprint(length(tm))\nlength(diag(TransM))')


# In[ ]:





# # Simulation of the ancestral history

# We can now simulate the ancestral history of 2 haploid samples under the coalescent with recombination using msprime.

# In[9]:


ts = msprime.sim_ancestry(
    samples=2,
    recombination_rate=0.1,
    sequence_length=100_000,
    random_seed = 237244,
    population_size = 1,
    ploidy = 1
)
# Visualise the simulated ancestral history.
# SVG(ts.draw_svg())


# We can now calculate the number of sites belonging to each tree inferred by msprime, together with the number of transitions between the trees: 

# In[10]:


# Create a matrix with 0 values
trans_mat = np.zeros((ts.num_nodes-2, ts.num_nodes-2))
first = True
# For each tree
for tree in ts.trees():
    # If not the first iteration
    if not first:
        # Count transition
        trans_mat[n_prev, tree.root-2] += 1
    # Update previous iteration counter
    n_prev = tree.root-2
    first = False
    # Add number of sites in that tree to diagonal
    trans_mat[tree.root-2, tree.root-2] += tree.span-1


# In[11]:


trans_mat[1:10, 1:10]


# We can also recover the height  each tree in the simulation:

# In[12]:


tree_heights = ts.tables.nodes[2:].asdict()['time'].tolist()
tree_heights[1:10]


# We now have all the ingredients to simplify the simulated transition matrix based on the cutpoints defined by the discretized times:

# In[13]:


get_ipython().run_cell_magic('R', '-i trans_mat -i tree_heights', '\n# Import the tre heights as vector\ntree_heights <- unlist(tree_heights)\n\n# Create vector for assigning each simulated tree to an interval\ncut_vec <- cut(tree_heights, breaks = c(0,tm), labels = FALSE, right = FALSE)\n\n# Create empty matrix with the right dimensions\ndiscrete_mat <- matrix(0, length(tm), length(tm))\n\n# For each interval\nfor (i in 1:length(tm)) {\n    # Save indices of trees belonging to interval i\n    vec_idx_i <- which(cut_vec == i)\n        # For each interval\n        for (j in 1:length(tm)) {\n            # Save indices of trees belonging to interval j\n            vec_idx_j <- which(cut_vec == j)\n            # If both indices are non-empty\n            if ((length(vec_idx_i) != 0) & (length(vec_idx_j) != 0)) {\n                # Save sum of the sub-matrix\n                discrete_mat[i, j] <- \n                    sum(trans_mat[vec_idx_i, vec_idx_j])   \n        }\n    }\n}\n\n\nprint(discrete_mat[1:10, 1:10])')


# In order to get the probabilities, we can divide by the total number per row. In theory, this should be the same as dividing it by the total number per column, but in the simulations there might be some differences due to the matrix not being completely symmetric. Just in case, we can compute the column sums after re-scaling by the row sums to check whether they are roughly the same:

# In[14]:


get_ipython().run_cell_magic('R', '', '\n# Divide by the row sums\ndiscrete_mat <- discrete_mat/rowSums(discrete_mat)\n\n# Check that  rowSums and colSums are roughly equal\nprint(rowSums(discrete_mat))\nprint(colSums(discrete_mat))')


# We can start comparing the theoretical and the simulated matrices. First, we can compare the diagonal values. 

# In[15]:


get_ipython().run_cell_magic('R', ' -w 700 -h 500', '\nggplot() +\n    geom_point(aes(tm, diag(discrete_mat), color = \'Empirical\')) +\n    geom_smooth(aes(tm, diag(discrete_mat), color = \'Empirical\')) +\n    geom_point(aes(tm, diag(TransM), color = \'Theoretical\')) +\n    geom_step(aes(tm, diag(TransM), color = \'Theoretical\'), direction="vh") +\n    theme_bw() +\n    theme(legend.position = \'bottom\') +\n    ylab(paste0(\'P(L = l, R = l)\')) +\n    xlab(paste0(\'Left tree height (l)\'))\n    ')


# Note that the x values for the points are the upper boundaries in tree height of each of the intervals. 

# We can also compute the probability given a certain interval. For example, we can plot the probabilities for when the left tree height is in the interval number 10:

# In[16]:


get_ipython().run_cell_magic('R', ' -w 700 -h 500', '\ninterval <- 10\n\nggplot() +\n    geom_point(aes(tm[-interval], discrete_mat[interval,][-interval], color = \'Empirical\')) +\n    geom_smooth(aes(tm[-interval], discrete_mat[interval,][-interval], color = \'Empirical\')) +\n    geom_step(aes(tm[-interval], TransM[interval,][-interval], color = \'Theoretical\'), direction="vh") +\n    geom_point(aes(c(0, tm[-interval]), c(0, TransM[interval,][-interval]), color = \'Theoretical\')) +\n    theme_bw() +\n    theme(legend.position = \'bottom\') +\n    ylab(paste0(\'P(L = l, R = r |\xa0d_l = \', interval, \')\')) +\n    xlab(paste0(\'Left tree height (l)\'))\n    ')


# Note that the x values for the points are the upper boundaries in tree height of each of the intervals. Also, the y value corresponding to the analyzed interval is also removed for simplicity. 

# We can also compare the matrices as a whole by plotting them side by side. The color scale is log-transformed to be able to see the values off the diagonal. Also, those values in the simulated matrix with a probability of 0 (because they were not found) are defaulted to the lowest color in the scale:

# In[17]:


get_ipython().run_cell_magic('R', ' -w 1000 -h 500', "\nmelt(discrete_mat) %>% \n    mutate(type = 'Simulation') %>%\n    bind_rows(\n        mutate(\n            melt(TransM),\n            type = 'Theoretical'\n        )\n    ) %>%\n    ggplot() +\n    geom_tile(aes(Var2, Var1, fill = value)) +\n    scale_y_reverse(breaks = 1:20, expand = c(0, 0)) +\n    scale_x_continuous(breaks = 1:20, expand = c(0, 0)) +\n    facet_wrap(~type) +\n    theme_bw() +\n    scale_fill_viridis_c(\n        trans = 'log', \n        na.value = hcl.colors(1)\n    )\n")


# # The mutational process

# After proving that the simulated values for the transition matrix match the theoretical expectation, we can now move on to analyzing the mutational process. We can calculate the probability for a mutation in a given discretized interval using the following:

# In[18]:


get_ipython().run_cell_magic('R', '', '\n\ncalc_p_mut <- function(tm, tht) {\n    tm_1 <- c(0, tm[1:(length(tm)-1)])\n    prob_mutation <- rep(NA, length(tm))\n    for (i in 1:length(tm)){\n        prob_mutation[i] <- 1-exp(-tht*tm_1[i])*(1-exp(-(1+tht)*(tm[i]-tm_1[i])))/(1-exp(-(tm[i]-tm_1[i])))/(1+tht)\n    }\n    prob_mutation\n}\n\nprob_mut <- calc_p_mut(tm, tht = 0.01*2)\n\nprob_mut')


# We are now ready to simulate the mutational process on top of the tree collection previously simulated using msprime:

# In[19]:


mutated_ts = msprime.sim_mutations(
    ts,
    rate = 0.01,
    model = 'infinite_alleles'
)

SVG(mutated_ts.first().draw_svg())


# We can now recover the number of mutations that happened in each of the trees:

# In[20]:


# Create a matrix with 0 values
mut_mat = np.zeros((mutated_ts.num_nodes-2, 2))
first = True
# For each tree
for tree in mutated_ts.trees():
    # Save the number of positions in the tree
    mut_mat[tree.root-2, 0] += tree.span
    # save the number of mutations
    mut_mat[tree.root-2, 1] += tree.num_mutations
    
mut_mat


# We can now calculate the number of mutations for each of the discretized intervals:

# In[21]:


get_ipython().run_cell_magic('R', '-i mut_mat', '\n(mutat_tib <- as_tibble(mut_mat) %>%\n    # add the interval associated to each tree\n    mutate(\n        bin = cut_vec\n    ) %>%\n    # group by each interval\n    group_by(bin) %>%\n    # calculate the number of sites \n    # and mutations per interval\n    summarize(\n        sites = sum(V1),\n        mutations = sum(V2),\n        p_obs = mutations/sites\n    ) %>%\n    # add the expected probability and height\n    mutate(\n        p_exp = prob_mut,\n        height = tm[bin]\n    ))')


# In[22]:


get_ipython().run_cell_magic('R', '-w 700 -h 500', '\nmutat_tib %>%    \n    pivot_longer(starts_with("p_")) %>%\n    ggplot() +\n        geom_point(aes(height, value, color = name)) +\n        geom_line(aes(height, value, color = name)) +\n        theme_minimal()')


# # Hidden Markov model

# We can think of the ancestral history with recombination and mutation as a hidden Markov model (HMM), where the discretized tree heights are the hidden states and the mutations are the emitted observations. We now have all the ingredients for building a hidden Markov model:
# 
# * The initial probabilities.
# * The transition probability matrix will be the transition matrix for the discretized times, which is dependant on rho.
# * The emission probability matrix will be the probability of observing (or not) a mutation, which is dependant on both rho and theta.

# The following code can be used to calculate the logarithmic forward likelihood given the initial probabilities, the transition probabilities, the emission probabilities and a vector of the observed mutations:

# In[23]:


get_ipython().run_cell_magic('R', '', '\nlogFrwdLikFct <- function(InitProb,TransProb,EndProb,EmisProb,ObsSeq){\n  # Number of observations\n  len <- length(ObsSeq)\n  # Number of hidden states\n  nHS <- nrow(TransProb)\n  # Define logForwardLik\n  logForwardLik <- matrix(0,nrow=len,ncol=nHS)\n  # Start condition\n  logForwardLik[1,] <- log(InitProb*EmisProb[,ObsSeq[1]])\n  # Determine logForwardLik by recursion \n  for(k in 2:len){\n    a <- max(logForwardLik[k-1,])\n    # Calculate the loglik of current iteration\n    logForwardLik[k,] <- \n      log(colSums(\n          (exp(logForwardLik[k-1,]-a)%*%TransProb)*EmisProb[,ObsSeq[k]]\n      ))+a\n  }  \n  a <- max(logForwardLik[len,])\n  logForwardLikVal <- log(sum(EndProb*exp(logForwardLik[len,]-a)))+a\n  return(logForwardLikVal)\n}')


# The above code has a little trick to improve numerical stability. For each iteration, it substracts the maximum loglik value of the previous iteration to the loglik values of the previous iteration before exponentiating them. After performing all the computations, we go back to the logarithm and we will add the value we substracted earlier.

# In order to run the forward algorithm, we should calculate the observed mutations first:

# In[24]:


mut_pos = []
for variant in mutated_ts.variants():
    mut_pos.append(int(variant.position))
mut_pos[1:10]


# We can create a binary vector, where a 1 means that there is a mutation in that position, and a 0 indicates a lack of mutation:

# In[25]:


mut_pos_bin = np.zeros(int(mutated_ts.sequence_length), dtype=np.int32)
for i in range(int(mutated_ts.sequence_length)):
    if i in mut_pos:
        mut_pos_bin[i] = 1
    else:
        mut_pos_bin[i] = 0
mut_pos_bin[1:10]


# ## HMM optimization over rho

# Let's now run the forward algorithm over a grid of rho values and save the loglik associated to them. In order to do so, we have to specify a fixed value of theta:

# In[26]:


get_ipython().run_cell_magic('R', '-i mut_pos_bin', '\ntht <- 0.01*2\nrho <- 0.1*2\n\nrho.v <- seq(from=0.10,to=0.30,by=0.01)\nlogL_rho.v <- rep(0,len=length(rho.v))\n\nfor (i in 1:length(rho.v)) {\n    initial <- rep(1/nInt, nInt)\n    transition <- FastTransMat(tm, rho.v[i])\n    emission <- matrix(c(1-calc_p_mut(tm, tht),calc_p_mut(tm, tht)),ncol=2,byrow=FALSE)\n    logL_rho.v[i] <- logFrwdLikFct(initial, transition, initial, emission, mut_pos_bin+1)\n    print(paste0("rho=",rho.v[i],", log likelihood=",logL_rho.v[i]))\n}\n')


# Instead of defining a grid of values ourselves, we can simply use the `optimize` function in R to calculate the value with the largest loglik:

# In[27]:


get_ipython().run_cell_magic('R', '', '\noptim_rho <- function(rho) {\n    initial <- rep(1/nInt, nInt)\n    transition <- FastTransMat(tm, rho)\n    emission <- matrix(c(1-calc_p_mut(tm, tht),calc_p_mut(tm, tht)),ncol=2,byrow=FALSE)\n    logFrwdLikFct(initial, transition, initial, emission, mut_pos_bin+1)\n}\n\nmax_lik_rho <- optimize(optim_rho, c(0.1, 0.3), maximum = TRUE)\n\nmax_lik_rho')


# When plotting these results together, we can see that the optimization worked quite well:

# In[28]:


get_ipython().run_cell_magic('R', '-w 700 -h 500', "\nggplot() +\n    geom_line(aes(rho.v, logL_rho.v)) +\n    geom_vline(aes(xintercept = rho, color = 'Real')) +\n    geom_vline(aes(xintercept = max_lik_rho$maximum, color = 'Optimization')) +\n    theme_minimal()")


# ## HMM optimization over theta

# Additionally, we can also optimize for theta by choosing a fixed value for rho:

# In[ ]:


get_ipython().run_cell_magic('R', '-i mut_pos_bin', '\ntht <- 0.01*2\nrho <- 0.1*2\n\ntht.v <- seq(from=0.01,to=0.03,by=0.001)\nlogL_tht.v <- rep(0,len=length(tht.v))\n\nfor (i in 1:length(tht.v)) {\n    initial <- rep(1/nInt, nInt)\n    transition <- FastTransMat(tm, rho)\n    emission <- matrix(c(1-calc_p_mut(tm, tht.v[i]),calc_p_mut(tm, tht.v[i])),ncol=2,byrow=FALSE)\n    logL_tht.v[i] <- logFrwdLikFct(initial, transition, initial, emission, mut_pos_bin+1)\n    print(paste0("tht=",tht.v[i],", log likelihood=",logL_tht.v[i]))\n}\n')


# In[ ]:


get_ipython().run_cell_magic('R', '', '\noptim_tht <- function(tht) {\n    initial <- rep(1/nInt, nInt)\n    transition <- FastTransMat(tm, rho)\n    emission <- matrix(c(1-calc_p_mut(tm, tht),calc_p_mut(tm, tht)),ncol=2,byrow=FALSE)\n    logFrwdLikFct(initial, transition, initial, emission, mut_pos_bin+1)\n}\n\nmax_lik_tht <- optimize(optim_tht, c(0.01, 0.03), maximum = TRUE)\n\nmax_lik_tht')


# In[ ]:


get_ipython().run_cell_magic('R', '-w 700 -h 500', "\nggplot() +\n    geom_line(aes(tht.v, logL_tht.v)) +\n    geom_vline(aes(xintercept = tht, color = 'Real')) +\n    geom_vline(aes(xintercept = max_lik_tht$maximum, color = 'Optimization')) +\n    theme_minimal()")


# ## HMM optimization over both rho and theta

# In[ ]:


get_ipython().run_cell_magic('R', '', "\nlibrary(optimParallel)\n\nfun_optim <- function(x) {\n    require(expm)\n    rho <- x[1]\n    tht <- x[2]\n    initial <- rep(1/nInt, nInt)\n    transition <- FastTransMat(tm, rho)\n    emission <- matrix(c(1-calc_p_mut(tm, tht),calc_p_mut(tm, tht)),ncol=2,byrow=FALSE)\n    logFrwdLikFct(initial, transition, initial, emission, mut_pos_bin+1)\n}\n\ncl <- makeCluster(20)     # set the number of processor cores\nsetDefaultCluster(cl=cl)  # set 'cl' as default cluster\nclusterExport(cl=cl, c('nInt', 'mut_pos_bin', 'tm', 'FastTransMat', 'calc_p_mut', 'logFrwdLikFct'))\noptim_both <- optimParallel(c(0.3, 0.05), \n                            fun_optim, \n                            lower = c(0.001, 0.001),\n                            upper = c(1, 1),\n                            parallel = list(loginfo = TRUE),\n                            control = list(ndeps=1e-4, fnscale = -1, pgtol = 0)\n                            )")


# In[ ]:


get_ipython().run_cell_magic('R', '', '\noptim_both')


# In[ ]:


get_ipython().run_cell_magic('R', '', '\noptim_both$loginfo %>%\n    as_tibble() %>%\n    mutate(\n        par1_end = lead(par1),\n        par2_end = lead(par2)) %>%\n    ggplot() +\n    geom_segment(aes(x = par1, y = par2, xend = par1_end, yend = par2_end), arrow = arrow()) +\n    geom_point(aes(par1, par2, color = -fn)) ')


# In[ ]:


get_ipython().system('jupyter nbconvert --to script 01_msprime_simulations.ipynb')

