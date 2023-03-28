
import numpy as np

def sample_xs_from_PTables(ptables_list):
    PTable = ptables_list[0]
    TotTable = ptables_list[1]
    CapTable = ptables_list[2]
    ScatTable = ptables_list[3]
    FisTable = ptables_list[4]

    rn = np.random.default_rng().uniform(low=0, high=1)
    xs_index = np.searchsorted(np.cumsum(PTable[0]), rn)
    Sig_t = TotTable[0,xs_index]
    Sig_g = CapTable[0, xs_index]
    Sig_s = ScatTable[0, xs_index]
    Sig_f = FisTable[0, xs_index]
 
    return Sig_t, Sig_g, Sig_s, Sig_f


def Eigen_function_0D_CE(E0, tally, material, rng, 
                                                URR_Erange = None,
                                                ptables_list = None,
                                                avg_URR = None):

    ### kills neutrons that exit lower energies
    if E0 < tally.Emin:
        return
    
    ### pull cross sections at E0
    if URR_Erange is None:
        Sig_t, Sig_f, Sig_g, Sig_s = material.get_macro_cross_sections(E0)
    elif ptables_list is None:
        assert avg_URR is not None
        if np.searchsorted(URR_Erange, E0) == 1:
            [Sig_t, Sig_g, Sig_s, Sig_f ]= avg_URR
        else:
            Sig_t, Sig_f, Sig_g, Sig_s = material.get_macro_cross_sections(E0)
    else:
        assert avg_URR is None
        if np.searchsorted(URR_Erange, E0) == 1:
            Sig_t, Sig_g, Sig_s, Sig_f = sample_xs_from_PTables(ptables_list)
        else:
            Sig_t, Sig_f, Sig_g, Sig_s = material.get_macro_cross_sections(E0)

    
    ### tally in energy
    tally.tally_energy(E0, Sig_t)

    ### sample reaction and do something
    reaction = rng.uniform(low=0.0, high=Sig_t)

    # if capture, return
    if reaction <= Sig_g:
        E_new = 0
        return
    
    # if fission, add neutrons and exit
    elif reaction <= Sig_g + Sig_f:
        nubar = material.nubar
        sample_nu = rng.uniform(low=np.floor(nubar), high=np.ceil(nubar))
        if sample_nu <= nubar:
            nu = np.ceil(nubar)
        else: # if sample_nu > nubar:
            nu = np.floor(nubar)  

        tally.first_moment += nu
        tally.second_moment += nu**2
        E_new = 0
        return 
    
    else:
        # %sample uniform from E0 to lower ebound
        E_new = rng.uniform(low=tally.Emin, high=E0)

        # % repeat energy transport function
        E_new = Eigen_function_0D_CE(E_new, tally, material, rng) ;  
        

    return E_new




def transport_loop_0D_CE(N, G, tally, mat, rng, URR_Erange, ptables_list, avg_URR):
        
    # option to set seed
    rng = np.random.default_rng()
    for g in range(int(G)):
        tally.reset_generation_tally()
        for iN in range(int(N)):
            #random fission energy
            E_start = rng.uniform(low=tally.Emin, high=tally.Emax) 
            # transport
            E_new = Eigen_function_0D_CE(E_start, tally, mat, rng,
                                                            URR_Erange = URR_Erange,
                                                            ptables_list = ptables_list,
                                                            avg_URR = avg_URR)

        # save generation tally
        tally.save_generation_tally(N)
        
    return tally




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