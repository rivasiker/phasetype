#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
from trails.optimizer import trans_emiss_calc
from trails.cutpoints import cutpoints_ABC, cutpoints_AB
import numpy as np
from trails.optimizer import forward_loglik, post_prob, viterbi
import pandas as pd
import time
import re
import msprime
get_ipython().run_line_magic('load_ext', 'rpy2.ipython')


# In[2]:


ILS = int(sys.argv[1])


# In[3]:


####################### Model parameters #######################

print('Specifying model')

N_AB = 50000
N_ABC = 30000
t_A = 160000
t_B = 190000
t_C = 110000
t_1 = max([t_A, t_B, t_C])
t_2 = -N_AB*np.log(3/2*ILS/100)
t_3 = t_1*5
r = 0.5e-8
mu = 2e-8
n_int_AB = 5
n_int_ABC = 7

t_out = t_1+t_2+t_3+2*N_ABC

N_ref = N_ABC

coal_ABC = N_ref/N_ABC
coal_AB = N_ref/N_AB
t_upper = t_3-cutpoints_ABC(n_int_ABC, coal_ABC)[-2]*N_ref
t_AB = t_2/N_ref

cut_AB = t_1+cutpoints_AB(n_int_AB, t_AB, coal_AB)*N_ref
cut_ABC = t_1+t_2+cutpoints_ABC(n_int_ABC, coal_ABC)*N_ref

(2/3)*(np.exp(-t_2/(N_AB)))


# In[4]:


print('Computing HMM')

transitions, emissions, starting, hidden_states, observed_states = trans_emiss_calc(
    t_A, t_B, t_C, t_2, t_upper, t_out,
    N_AB, N_ABC,
    r, mu, mu, mu, mu, mu, mu, n_int_AB, n_int_ABC)

dct_hid = {v: k for k, v in hidden_states.items()}
dct = {v: k for k, v in observed_states.items()}


# In[5]:


####################### Add demography #######################

print('Adding demography')

n_sites = 200_000
seed = 10

demography = msprime.Demography()
demography.add_population(name="A", initial_size=N_AB, default_sampling_time=t_1-t_A)
demography.add_population(name="B", initial_size=N_AB, default_sampling_time=t_1-t_B)
demography.add_population(name="C", initial_size=N_AB, default_sampling_time=t_1-t_C)
demography.add_population(name="D", initial_size=N_AB, default_sampling_time=t_1-t_1)
demography.add_population(name="AB", initial_size=N_AB)
demography.add_population(name="ABC", initial_size=N_ABC)
demography.add_population(name="ABCD", initial_size=N_ABC)
demography.add_population_split(time=t_1, derived=["A", "B"], ancestral="AB")
demography.add_population_split(time=t_1+t_2, derived=["AB", "C"], ancestral="ABC")
demography.add_population_split(time=t_1+t_2+t_3, derived=["ABC", "D"], ancestral="ABCD")

ts = msprime.sim_ancestry(
    {"A": 1, "B": 1, "C": 1,
     "D": 1
    },
    demography=demography,
    recombination_rate=r*2,
    sequence_length=n_sites,
    ploidy=1,
    random_seed=seed
)


# In[6]:


#### Add mutations

print('Adding mutations')

mutated_ts = msprime.sim_mutations(ts, rate=mu*2, random_seed=seed)

nochange_lst = [dct['AAAA'], dct['CCCC'], dct['TTTT'], dct['GGGG']]
np.random.seed(seed) ; sim_genome = np.random.choice(nochange_lst, n_sites)

mut_lst = []
mut_loc = []
for variant in mutated_ts.variants():
    mut_loc.append(variant.site.position)
    mut_lst.append(''.join([variant.alleles[i] for i in variant.genotypes]))

for i in range(len(mut_loc)):
    sim_genome[int(mut_loc[i])] = dct[mut_lst[i]]


# In[7]:


# loglik = forward_loglik(transitions, emissions, starting, sim_genome)


# In[8]:


post = post_prob(transitions, emissions, starting, sim_genome)


# In[9]:


# vit = viterbi(transitions, emissions, starting, sim_genome)


# In[10]:


hidden_matrix = np.random.randint(max([n_int_AB, n_int_ABC]), size=(len(dct_hid), 4))
hidden_matrix[:,0] = list(range(len(dct_hid)))
hidden_matrix[:,1] = [i[0] for i in dct_hid.keys()]
hidden_matrix[:,2] = [i[1] for i in dct_hid.keys()]
hidden_matrix[:,3] = [i[2] for i in dct_hid.keys()]


# In[11]:


print('Computing cuts')

left_lst = []
right_lst = []
tree_state = []
t_AB_vec = []
t_ABC_vec = []
for t in ts.trees():
    # Append start coordinate
    left_lst.append(t.interval.left)
    # Append end coordinate
    right_lst.append(t.interval.right-1)
    # Get all non-zero coalescent times
    ntimes = [ts.nodes()[n].time for n in t.nodes() if ts.nodes()[n].time not in [0, t_1-t_A, t_1-t_B, t_1-t_C]]
    ntimes = sorted(ntimes)
    # Get time of the first event
    mint = ntimes[0]
    mint2 = ntimes[1]
    # Find topology
    find_re = re.findall("n\d,n\d", t.as_newick(include_branch_lengths=False))[0]
    # Sort species within topology
    find_re = sorted(find_re.split(','))
    # If V0 or V1
    if find_re == ['n0', 'n1']:
        # If the time of the first coalescent is larger than the deepest speciation event
        if mint>=(t_1+t_2):
            state = (1, (mint>cut_ABC).sum()-1, (mint2>cut_ABC).sum()-1)
            # Append V1 state
        else:
            state = (0, (mint>cut_AB).sum()-1, (mint2>cut_ABC).sum()-1)
            # Append V0 state
    # If V2
    elif find_re == ['n0', 'n2']:
        state = (2, (mint>cut_ABC).sum()-1, (mint2>cut_ABC).sum()-1)
    # If V3
    elif find_re == ['n1', 'n2']:
        state = (3, (mint>cut_ABC).sum()-1, (mint2>cut_ABC).sum()-1)
    else:
        state = (4, (mint>cut_ABC).sum()-1, (mint2>cut_ABC).sum()-1)
    tree_state.append(state)
    t_AB_vec.append(mint)
    t_ABC_vec.append(mint2)


# In[12]:


tree_matrix = np.random.randint(max(left_lst), size=(len(left_lst), 3))
tree_matrix[:,0] = left_lst
tree_matrix[:,1] = right_lst
tree_matrix[:,2] = [dct_hid[i] for i in tree_state]

print('Running R')


# In[13]:


get_ipython().run_cell_magic('R', '-i post -i hidden_matrix -i tree_matrix', "\nlibrary(tidyverse)\n\nhid_tab <- as_tibble(hidden_matrix) %>%\n    rename(name = V1, topology = V2, int_1 = V3, int_2 = V4)\n    \n# write_csv(hid_tab, 'hid_tab.csv')\n    \ntree_tab <- as_tibble(tree_matrix) %>%\n    rename(start = V1, end = V2, name = V3) %>%\n    mutate(\n        gr = ifelse(lag(name) != name, 1, 0) %>% coalesce(0),\n        gr = cumsum(gr) + 1\n    ) %>% \n    group_by(gr, name) %>%\n    summarize(start = min(start), end = max(end)) %>%\n    left_join(hid_tab, by = 'name')\n    \n# write_csv(tree_tab, 'tree_tab.csv')\n\npost_tab <- as_tibble(post) %>%\n    mutate(pos = 0:(n()-1)) %>%\n    pivot_longer(-pos) %>%\n    mutate(name = as.integer(str_remove_all(name, 'V'))-1) %>%\n    left_join(hid_tab, by = 'name')\n     \n# write_csv(post_tab, 'post_tab.csv')")


# In[15]:


get_ipython().run_cell_magic('R', '-w 2000 -h 700 -r 150 -i n_int_AB', '\nd1 <- post_tab %>%\n    mutate(is_V0 = topology == 0) %>%\n    group_by(pos, is_V0, int_1) %>%\n    summarize(prob = sum(value))\nd2 <- post_tab %>%\n    group_by(pos, int_2) %>%\n    summarize(prob = sum(value))\nd3 <- post_tab %>%\n    group_by(pos, topology) %>%\n    summarize(prob = sum(value)) \n    \nmin_prob <- min(c(d1$prob, d2$prob, d3$prob))\nmax_prob <- max(c(d1$prob, d2$prob, d3$prob))')


# In[16]:


get_ipython().run_cell_magic('R', '-w 2000 -h 700 -r 150', "\np1 <- d1 %>%\n    ggplot() +\n    geom_tile(aes(pos, int_1+(!is_V0)*(n_int_AB+0.1), fill = prob, color = prob)) +\n    geom_segment(aes(x = start, xend = end, y = int_1+(!(topology == 0))*(n_int_AB+0.1), yend = int_1+(!(topology == 0))*(n_int_AB+0.1)), \n                 color = 'red',\n                 data = tree_tab) +\n    scale_fill_viridis_c(limits=c(0,1)) +\n    scale_color_viridis_c(limits=c(0,1)) +\n    scale_x_continuous(expand = c(0, 0)) +\n    scale_y_continuous(expand = c(0, 0)) +\n    labs(y = 'First coalescent', x = 'Position')\n\np2 <- d2 %>%\n    ggplot() +\n    geom_tile(aes(pos, int_2, fill = prob, color = prob)) +\n    geom_segment(aes(x = start, xend = end, y = int_2, yend = int_2), \n                 color = 'red',\n                 data = tree_tab) +\n    scale_fill_viridis_c(limits=c(0,1)) +\n    scale_color_viridis_c(limits=c(0,1)) +\n    scale_x_continuous(expand = c(0, 0)) +\n    scale_y_continuous(expand = c(0, 0)) +\n    labs(y = 'Second coalescent') +\n    theme(\n        axis.title.x = element_blank(),\n        axis.text.x = element_blank(),\n        axis.ticks.x=element_blank()\n    )\n    \np3 <- d3 %>%\n    ggplot() +\n    geom_tile(aes(pos, topology, fill = prob, color = prob)) +\n    geom_segment(aes(x = start, xend = end, y = topology, yend = topology), \n                 color = 'red',\n                 data = tree_tab) +\n    scale_fill_viridis_c(limits=c(0,1)) +\n    scale_color_viridis_c(limits=c(0,1)) +\n    scale_x_continuous(expand = c(0, 0)) +\n    scale_y_continuous(expand = c(0, 0)) +\n    labs(y = 'Topology') +\n    theme(\n        axis.title.x = element_blank(),\n        axis.text.x = element_blank(),\n        axis.ticks.x=element_blank()\n    )")


# In[18]:


get_ipython().run_cell_magic('R', '-w 2000 -h 1000 -r 150 -i ILS', "\nlibrary(patchwork)\n\np3/p2/p1+plot_annotation(\n  title = paste0('Posterior decoding, ILS = ', ILS, '%')\n)\n\nggsave(paste0('plot_posterior_decoding_', round(ILS), '.png'), width = 14, height = 7)")


# In[ ]:




