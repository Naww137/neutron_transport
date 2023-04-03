import numpy as np
import pandas as pd


def sample_spacing(D):
    spacing = np.sqrt(-4/np.pi*np.log( np.random.default_rng().uniform(low=0.0,high=1.0,size=2) ))
    return spacing*D

# def calculate_D3_statistic(E_levels):
#     Ef = E_levels[-1]
#     Ei = E_levels[0]
#     N = []
#     Erange = np.linspace(Ei, Ef, 100)
#     for El in Erange:
#         test = 0
#         N.append(np.count_nonzero(np.array(E_levels)[E_levels<El]))
#     def linear(vars):
#         return vars[0]*Erange+vars[1]
#     def fit(vars):
#         return np.sum( (linear(vars) - np.array(N))**2)
#     out = minimize(fit, x0=(1e3,0))
#     D3 = out.fun/2/len(E_levels)
#     D3_theo = 1/np.pi**2 * np.log(len(E_levels)-0.06871) 
#     print(f'The final D3 statistic and theoretical value are {D3} and {D3_theo}')
#     return

def get_resonance_ladder(E_ref, resonance_pairs, D_avg, Gg_avg, Gn_avg, Gf_avg):

    pairs = 0
    rn = np.random.default_rng().uniform(low=0.0,high=1.0)

    El_pos = [E_ref+rn*D_avg]
    El_neg = [E_ref+(rn-1)*D_avg]

    while pairs < resonance_pairs:
        spacing = sample_spacing(D_avg)
        El_neg.insert(0, El_neg[0]-spacing[0])
        El_pos.append(El_pos[-1]+spacing[1])
        E_levels = El_neg + El_pos
        # calculate D3 statistic to determine when to stop
        pairs += 1

    # sample reaction widths
    # Gg = np.random.chisquare(1000, size=len(E_levels))
    Gg = np.array([Gg_avg]*len(E_levels))
    Gn =  np.array([Gn_avg]*len(E_levels))*np.random.chisquare(1, size=len(E_levels))
    Gf =  np.array([Gf_avg]*len(E_levels))*np.random.chisquare(3, size=len(E_levels))

    E_levels = np.array(E_levels)
    # calculate_D3_statistic(E_levels)
    ladder = pd.DataFrame({'E'    :    E_levels,
                            'Gn'    :   Gn,
                            'Gg'    :   Gg,
                            'Gf'    :   Gf})
    
    return ladder

def calculate_xs_history(energy_grid, E_ref, isotope_or_mat):
    rxns = isotope_or_mat.get_macro_cross_sections(energy_grid)
    rxns_Eref = [np.unique(each[energy_grid==E_ref]) for each in rxns] 
    # np.unique([energy_grid==E_ref]) tot, cap, scat, fis
    return rxns_Eref, rxns