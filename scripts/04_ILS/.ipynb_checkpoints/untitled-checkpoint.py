import sys
import numpy as np
import pandas as pd
from scipy.linalg import expm
from scipy.stats import truncexpon
from scipy.stats import expon
from scipy.special import comb
import ast
import multiprocess as mp

sys.path.append('./functions')
from cutpoints_ABC import cutpoints_ABC
from get_ABC_inf_bis import get_ABC_inf_bis
from vanloan_3 import vanloan_3
from get_times import get_times
from get_tab_AB import get_tab_AB
from get_ordered import get_ordered
from vanloan_2 import vanloan_2
from cutpoints_AB import cutpoints_AB
from instant_mat import instant_mat
from vanloan_1 import vanloan_1
from get_ABC import get_ABC
from get_tab_ABC import get_tab_ABC
from get_joint_prob_mat import get_joint_prob_mat
from combine_states import combine_states
from load_trans_mat import load_trans_mat
from trans_mat_num import trans_mat_num

# Load variables
t_A = float(sys.argv[1])
t_B = float(sys.argv[2])
t_AB = float(sys.argv[3])
t_C = float(sys.argv[4])
rho_A = float(sys.argv[5])
rho_B = float(sys.argv[6])
rho_AB = float(sys.argv[7])
rho_C = float(sys.argv[8])
rho_ABC = float(sys.argv[9])
coal_A = float(sys.argv[10])
coal_B = float(sys.argv[11])
coal_AB = float(sys.argv[12])
coal_C = float(sys.argv[13])
coal_ABC = float(sys.argv[14])
n_int_AB = int(sys.argv[15])
n_int_ABC = int(sys.argv[16])

# Run function
array_np = get_joint_prob_mat(
    t_A,    t_B,    t_AB,    t_C, 
    rho_A,  rho_B,  rho_AB,  rho_C,  rho_ABC,
    coal_A, coal_B, coal_AB, coal_C, coal_ABC,
    n_int_AB, n_int_ABC)

# 0.1, 0.1, 1, 2, 2,   1,   3, 1, 1, 1,   0.5, 1, 1, 1, 1, 1

print(array_np[:,2].sum())