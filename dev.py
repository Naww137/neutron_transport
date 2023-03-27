# %%
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

from transport import functions_for_scattering_theory as xs


# %%
# %matplotlib widget

# %% [markdown]
# ### Match Benchmark
# Plot and print out values
# 
# 
# Incorrect statement:
# "Since the value of ν is 0 for the capture isotope, the total fission cross section is given only by the fission isotope."

# %%
resonance_ladder_pu9 = pd.DataFrame({'E'    :   [2.956242e-1],
                                'Gn'    :   [7.947046e-5],
                                'Gg'    :   [3.982423e-2],
                                'Gf'    :   [5.619673e-2]})

pair_constants_pu9 = {
'ac':9.41e-4,             # 1e-12cm
'p0':0.002196807122623 * 1/2 ,   # 1 / 1e-12cm - sqrt(eV)
'gj':3/4,
}

resonance_ladder_u8 = pd.DataFrame({'E'    :    [6.674280e10],
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



# %%

benchmark_e = np.array([0.00001, 0.01000, 0.29562, 6.67425, 100.000, 2.0e7])
# benchmark_e = np.array([6.67425])
energy = np.logspace(-5,7, 1000)
# energy = benchmark_e


pu9_Sig_t, pu9_Sig_f, pu9_Sig_g, pu9_Sig_s = pu9.get_macro_cross_sections(energy)
u8_Sig_t, u8_Sig_f, u8_Sig_g, u8_Sig_s = u8.get_macro_cross_sections(energy)

m_Sig_t, m_Sig_f, m_Sig_g, m_Sig_s = mat.get_macro_cross_sections(energy)


if np.array_equal(energy,benchmark_e):
    print(f'fission:\n {m_Sig_f}')
    print(f'scattering:\n {m_Sig_s}')
    print(f'total:\n {m_Sig_t}')
else:
    figure()
    # plot(energy, m_Sig_f, label='fission')
    # # plot(energy, mat_gam, label='gamma')
    # plot(energy, m_Sig_s, label='scattering')
    # plot(energy, m_Sig_t, label='total')

    plot(energy, u8_Sig_t, label='U8 total')
    plot(energy, u8_Sig_f, label='U8 f')

    xscale('log'); yscale('log')
    legend()

# %% [markdown]
# ## Sample a URR Ladder
# For Pu-239
# For U238, the URR occurs from 20 keV to just above 149 keV. 
# Let's start simple, get rid of the Pu-239 isotope and do the following:
# 1. Take a reference point $E_{ref}$ = 100 keV
# 2. Generate resonance-pair sequences out to either side of $E_{ref}$ and reconstruct the Fission, Capture, and Elastic scattering cross sections. (Make a nice plot!)
# 3. Use the Dyson-Mehta $\Delta_3$ statistic to terminate the end of the simulated resonance ladder

# %%
u8_Sig_t, u8_Sig_f, u8_Sig_g, u8_Sig_s = u8.get_macro_cross_sections(energy)

D_avg = 4
Gn_avg = 1.492300e-3
Gg_avg = 2.271100e-2
Gf_avg = 9.880000e-9

Ebins = np.array([100])*1e2
vEbins = np.array([20,149])*1e2

def sample_spacing(D):
    spacing = np.sqrt(-4/np.pi*np.log( np.random.default_rng().uniform(low=0.0,high=1.0,size=2) ))
    return spacing*D

# for E ref in Ebins
E_ref = Ebins[0]
steps = 0

rn = np.random.default_rng().uniform(low=0.0,high=1.0)

El_pos = [E_ref+rn*D_avg]
El_neg = [E_ref+(rn-1)*D_avg]

while steps < 500:

    spacing = sample_spacing(D_avg)
    El_neg.insert(0, El_neg[0]-spacing[0])
    El_pos.append(El_pos[-1]+spacing[1])
    E_levels = El_neg + El_pos
    # print(E_levels)

    # calculate D3 statistic to determine when to stop
    steps += 1

# sample reaction widths
# Gg = np.random.chisquare(1000, size=len(E_levels))
Gg = np.array([Gg_avg]*len(E_levels))
Gn = np.random.chisquare(1, size=len(E_levels))
Gf = np.random.chisquare(3, size=len(E_levels))

URR_energy = np.linspace(20e2, 149e2, 10000)
E_levels = np.array(E_levels)

u8_URR = pd.DataFrame({'E'    :    E_levels,
                        'Gn'    :   Gn,
                        'Gg'    :   Gg,
                        'Gf'    :   Gf})

u8_URR = material_definitions.isotope(pair_constants_u8, u8_URR, Nu8, 0)
u8_Sig_t, u8_Sig_f, u8_Sig_g, u8_Sig_s = u8_URR.get_macro_cross_sections(URR_energy)

figure()
plot(energy, m_Sig_t, label='total')
plot(URR_energy, u8_Sig_g)
plot(URR_energy, u8_Sig_s)
plot(URR_energy, u8_Sig_f)
plot(URR_energy, u8_Sig_t, 'b')
xscale('log')
yscale('log')

# print(URR_energy[u8_Sig_t<u8_Sig_f])
# u8_Sig_t, u8_Sig_f, u8_Sig_g, u8_Sig_s = u8_URR.get_macro_cross_sections(URR_energy[u8_Sig_t<u8_Sig_f])

#%%


# %% [markdown]
# #### Need to figure out the dyson mehta statistic for when to quit sampling resonances!

# %%
## calculate D3 statistic
# Ef = E_levels[-1]
# Ei = E_levels[0]
# N = []
# Erange = np.linspace(Ei, Ef, 100)
# for El in Erange:
#     test = 0
#     N.append(np.count_nonzero(np.array(E_levels)[E_levels<El]))

# figure()
# plot(Erange,N)
# # xscale('log')

# def linear(vars):
#     return vars[0]*Erange+vars[1]
# def fit(vars):
#     return np.sum( (linear(vars) - np.array(N))**2)

# from scipy.optimize import minimize
# out = minimize(fit, x0=(1e3,0))

# print(out.fun/2/len(E_levels))
# print(1/np.pi**2 * np.log(len(E_levels)-0.06871) )

# plot(Erange, linear(out.x))







# iL = np.arange(1, len(E_levels))
# N_E = np.sum(iL**2 * (np.array(E_levels[1:]) - np.array(E_levels[:-1])))

# g1 = np.sum( iL * (np.array(E_levels[1:])**2 - np.array(E_levels[:-1])**2) /2 )
# g2 = np.sum( iL * (np.array(E_levels[1:]) - np.array(E_levels[:-1])) )

# a1 = (Ef**3-Ei**3)/3
# a2 = (Ef**2-Ei**2)/2
# b1 = a2
# b2 = Ef-Ei

# a = (g1-g2)*(b1/b2) / (a1-a2)*(b1/b2)
# b = g2/b2 - a2/b2 * a

# D3 = 1/(Ef-Ei)*N_E - g1*a - g2*b

# from scipy.optimize import minimize

# E_levels = np.array(E_levels)
# iL = np.arange(1, len(E_levels))

# # def Dfunc(Ef, Ei, E_levels, iL, a,b):
# #     D = 1/len(E_levels) * np.sum( (iL*(E_levels[1:]-E_levels[:-1]) - a*(Ef**2-Ei**2)/2 - b)**2 )
# #     return D

# def Dfunc(X, XB, A:float, B:float):
#     N  = len(X)
#     H  = np.arange(N+1)
#     Y  = A*X+B
#     PB = (A*XB[0]+B, A*XB[1]+B-N)
#     P1 = Y-H[:-1]
#     P2 = Y-H[1:]
#     return (np.sum(P1**3 - P2**3) + (PB[1]**3 - PB[0]**3))/(3*A*(XB[1]-XB[0]))

# def func(indvars):  
#     return Dfunc(E_levels, vEbins, *indvars)

# sol  = minimize(func, x0=(1,1))

# print(sol)
# print( 1/np.pi**2 * np.log(len(E_levels)-0.06871) )


# %%



# %%
test_list = [4, 5, 6, 3, 9]
insert_list = [2, 3]
 
# initializing position
pos = 2
 
# printing original list
print ("The original list is : " + str(test_list))
 
# printing insert list
print ("The list to be inserted is : " + str(insert_list))
 
# using list slicing
# to insert one list in another
test_list[pos:pos] = insert_list
 
# printing result
print ("The list after insertion is : " + str(test_list))

# %%



# %%


# %%


# %%


# %%
#### problem 4
N = 1.5e2
G=20
iE = 500
print(f'Particle Histories: {N}')
print(f'Energy Bins: {iE}')

Emin = 1e-5
Emax = 2e7 
iEbins = 200
tally = tallies.tallies(Emin, Emax, iEbins)

# option to set seed
rng = np.random.default_rng()


# %%
#### Could try to run as a vector - later when I want to speed this up
# E0 = rng.uniform(low=tally.Emin, high=tally.Emax, size=(1,2))
# E0 < np.full_like(E0, tally.Emin)
# E_new = Eigen_function_0D_CE(E0, tally, mat, rng)


# %% [markdown]
# #### Run transport

# %%
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

# %%
figure()
errorbar(tally.Ebins, collision_based_scalar_flux[0], yerr=np.sqrt(collision_based_scalar_flux[1]), 
                                                            fmt='.', color='k', ms=3, ecolor='r', capsize=2)

xscale('log')
yscale('log')


# %%
tally.first_moment/N

# %%


# %%


# %%



