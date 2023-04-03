#%%

import numpy as np
import numpy as np
import pandas as pd
from matplotlib.pyplot import *
import numba as nb
import sys
import os

sys.path.insert(0, os.path.realpath('./transport/'))

from transport import tallies
from transport import material_definitions
from transport import functions_for_transport


i = int(sys.argv[1])
#%% original benchmark xs

resonance_ladder_pu9 = pd.DataFrame({'E'    :   [2.956242e-1], #, 2.5e3, 3e3],
                                'Gn'    :   [7.947046e-5], #, 7.947046e-5, 7.947046e-5],
                                'Gg'    :   [3.982423e-2], #, 3.982423e-2, 3.982423e-2],
                                'Gf'    :   [5.619673e-2]}) #, 5.619673e-2, 5.619673e-2]})

pair_constants_pu9 = {
'ac':0,             # 1e-12cm
'p0':0.002196807122623 * 1/2 ,   # 1 / 1e-12cm - sqrt(eV)
'gj':3/4,
}

resonance_ladder_u8 = pd.DataFrame({'E'    :    [6.674280], # 2.5e3, 3e3],
                                    'Gn'    :   [1.492300e-3], # 1.492300e0, 1.492300e0],
                                    'Gg'    :   [2.271100e-2], # 2.271100e0, 2.271100e0],
                                    'Gf'    :   [9.880000e-9]}) # 9.880000e-9, 9.880000e-9]})

pair_constants_u8 = {
'ac':0,   # 1e-12cm
'p0':0.002196807122623 * 1/2 ,   # 1 / 1e-12cm - sqrt(eV)
'gj':1,
}

Npu9 = 1
Nu8 = 0.124954
Nscat = 0.008340505

pu9 = material_definitions.isotope(pair_constants_pu9, resonance_ladder_pu9, Npu9, 2.88)
u8 = material_definitions.isotope(pair_constants_u8, resonance_ladder_u8, Nu8, 0)

mat_no_URR = material_definitions.material([pu9,u8], constant_scattering=20*Nscat)

#%% Load URR xs at E_ref history

xs_histories = np.genfromtxt('./xs_histories.csv', delimiter=',')


#%% Run transport for history

# # test
N = 1.5e2
G = 200
iE = 200
print(f'Particle Histories: {N}')
print(f'Energy Bins: {iE}')

Emin = 1e-5
Emax = 2e7 
iEbins = 200
tally = tallies.tallies(Emin, Emax, iEbins)

URR_Erange = np.array([20,149])*1e3
URR_single_history_xs = [xs_histories[i,0], xs_histories[i,1], xs_histories[i,2], xs_histories[i,3]]
rng = np.random.default_rng()

tally_avg = functions_for_transport.transport_loop_0D_CE(N, G, tally, mat_no_URR, rng, URR_Erange, None, URR_single_history_xs)
final_k_estimate, final_estimator_variance, collision_based_scalar_flux = tally.final_analysis() 

print(final_k_estimate, final_estimator_variance)




# %%

with open(f'transport_result_{i}.csv', 'w') as f:
    f.write(f'{final_k_estimate}, {final_estimator_variance}\n')
    for each in collision_based_scalar_flux.flatten():
        f.write(f'{each},')


# %%


# %%
