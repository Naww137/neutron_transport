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

from transport.functions_for_scattering_theory import sigma_s, sigma_g, sigma_f

#%%

#### problem 4
N = 1.5e3
G=3
iE = 200
print(f'Particle Histories: {N}')
print(f'Energy Bins: {iE}')

Emin = 1e-5
Emax = 2e7 
iEbins = 200
tally = tallies.tallies(Emin, Emax, iEbins)


# option to set seed
rng = np.random.default_rng()

#%% define isotopes and material

resonance_ladder_pu9 = pd.DataFrame({'E'    :   [2.956242e-1],
                                'Gn'    :   [7.947046e-5],
                                'Gg'    :   [3.982423e-2],
                                'Gf'    :   [5.619673e-2]})

pair_constants_pu9 = {
'ac':9.41e-4,             # 1e-12cm
'p0':0.002196807122623 * 1/2 ,   # 1 / 1e-12cm - sqrt(eV)
'gj':3/4,
}

resonance_ladder_u8 = pd.DataFrame({'E'    :    [6.674280],
                                    'Gn'    :   [1.492300e-3],
                                    'Gg'    :   [2.271100e-2],
                                    'Gf'    :   [9.880000e-9]})

pair_constants_u8 = {
'ac':9.48e-4,             # 1e-12cm
'p0':0.002196807122623 * 1/2 ,   # 1 / 1e-12cm - sqrt(eV)
'gj':1,
}

Npu9 = 1
Nu8 = 0.124954
Nscat = 0.008340505

pu9 = material_definitions.isotope(pair_constants_pu9, resonance_ladder_pu9, Npu9, 2.88)
u8 = material_definitions.isotope(pair_constants_u8, resonance_ladder_u8, Nu8, 0)

mat = material_definitions.material([pu9,u8], constant_scattering=20*Nscat)

        
#%%


for g in range(int(G)):

    tally.reset_generation_tally()
    
    for iN in range(int(N)):
        #random fission energy
        E_start = rng.uniform(low=tally.Emin, high=tally.Emax) 
        # transport
        E_new = functions_for_transport.Eigen_function_0D_CE(E_start, tally, mat, rng)
    
    # save generation tally
    tally.save_generation_tally(N)
    
final_k_estimate, final_estimator_variance, collision_based_scalar_flux = tally.final_analysis() 

print(final_k_estimate, final_estimator_variance)


# %%
figure()
plot(tally.Ebins, collision_based_scalar_flux[0])

xscale('log')
yscale('log')
show()
close()


# %%
