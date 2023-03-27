
import numpy as np





def Eigen_function_0D(E0, tally, material, rng):

    Sig_t, Sig_f, Sig_g, Sig_s = material.get_macro_cross_sections_groupwise()

    # define new E value
    E_new = E0 

    reaction = rng.uniform(low=0.0, high=Sig_t)
    
    # if capture, exit
    if reaction <= Sig_g: 
        return
    # if fission, add neutrons and exit
    elif reaction <= Sig_g + Sig_f:

        sample_nu = rng.uniform(low=0.0, high=1.0)
        if sample_nu <=0.6:
            nu = 2
        elif sample_nu > 0.6:
            nu = 3
            
        tally.first_moment += nu
        tally.second_moment += nu**2
        return
    # if scattering, run the transport again starting from collision location (in E)
    else:
        E_new = Eigen_function_0D(E_new, tally, material, rng)

    return E_new