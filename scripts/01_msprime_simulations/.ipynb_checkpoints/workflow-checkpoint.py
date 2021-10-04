from gwf import Workflow
import os
import shutil

gwf = Workflow()


rho_v = [0.0001, 0.005, 0.01, 0.05, 0.1, 0.5]
tht_v = [0.0001, 0.005, 0.01, 0.05, 0.1, 0.5]



# For each species
for rho in rho_v:
    for tht in tht_v:
        for ide in range(1, 6, 1): 
            gwf.target('sim_{}_{}_{}'.format(ide, rho, tht),
                        inputs=[], 
                        outputs=['../../steps/01_msprime_simulations/sim{}_rho{}_tht{}.csv'.format(ide, rho, tht)],
                        cores=6,
                        memory='4g',
                        walltime= '01:00:00',
                        account='Primategenomes') << """
                ipython 02_msprime_simulations.py {} {} {}
                """.format(ide, rho, tht)
 