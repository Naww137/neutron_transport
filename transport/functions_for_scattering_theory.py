

import numpy as np
import pandas as pd



def unpack_resonance_ladder(resonance_ladder):
    Elam = np.array(resonance_ladder.E)
    Gn = np.array(resonance_ladder.Gn)
    Gg = np.array(resonance_ladder.Gg)
    Gf = np.array(resonance_ladder.Gf)
    return Elam, Gn, Gg, Gf

def unpack_pair_constants(pair_constants):
    ac = pair_constants['ac']
    gj = pair_constants['gj']
    p0 = pair_constants['p0']
    return ac, gj, p0 

def Gn_vector(E, Elam, Gn):
    return Gn*np.sqrt(E/Elam)

def d(E, Elam, Gt_vec):
    return (E-Elam)**2 + (Gt_vec/2)**2

def sigma_g(E, pair_constants, resonance_ladder):

    ac, gj, p0 = unpack_pair_constants(pair_constants)
    k = p0*np.sqrt(E)

    xs = 0
    for index, row in resonance_ladder.iterrows():
        Elam, Gn, Gg, Gf = unpack_resonance_ladder(row)
        Gn_vec = Gn_vector(E, Elam, Gn)
        Gt_vec = Gn_vec+Gg+Gf

        xs +=  Gn_vec * Gg / d(E, Elam, Gt_vec)
    xs = gj*np.pi/k**2 * xs
    return xs

def sigma_f(E, pair_constants, resonance_ladder):

    ac, gj, p0 = unpack_pair_constants(pair_constants)
    k = p0*np.sqrt(E)

    xs = 0
    for index, row in resonance_ladder.iterrows():
        Elam, Gn, Gg, Gf = unpack_resonance_ladder(row)
        
        Gn_vec = Gn_vector(E, Elam, Gn)
        Gt_vec = Gn_vec+Gg+Gf

        xs += Gn_vec * Gf / d(E, Elam, Gt_vec)
    xs = gj*np.pi/k**2 * xs
    return xs

def sigma_s(E, pair_constants, resonance_ladder):

    ac, gj, p0 = unpack_pair_constants(pair_constants)
    k = p0*np.sqrt(E)

    xs = 0
    for index, row in resonance_ladder.iterrows():
        Elam, Gn, Gg, Gf = unpack_resonance_ladder(row)

        Gn_vec = Gn_vector(E, Elam, Gn)
        Gt_vec = Gn_vec+Gg+Gf

        xs += gj*np.pi/d(E, Elam, Gt_vec) * ( Gn_vec**2/k**2 + 4*ac*(E-Elam)*Gn_vec/k - 2*ac**2*Gn_vec*Gt_vec )
    xs = 4*np.pi*ac**2 + xs
    return xs