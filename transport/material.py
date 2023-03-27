
import numpy as np
import functions_for_scattering_theory as xs


class material:

    def __init__(self, pair_constants, resonance_ladder, number_density):

        self.pair_constants = pair_constants
        self.resonance_ladder = resonance_ladder
        self.number_density = number_density

        return
    

    def get_micro_cross_sections(self, energy):

        sig_g = xs.sigma_g(energy, self.pair_constants, self.resonance_ladder)
        sig_s = xs.sigma_s(energy, self.pair_constants, self.resonance_ladder)
        sig_f = xs.sigma_f(energy, self.pair_constants, self.resonance_ladder)
        sig_t = sig_g+sig_f+sig_s

        return sig_t, sig_f, sig_g, sig_s

    def get_macro_cross_sections(self, energy):

        sig_t, sig_f, sig_g, sig_s = self.get_micro_cross_sections(energy)
        Sig_t, Sig_f, Sig_g, Sig_s = self.number_density*np.array([sig_t, sig_f, sig_g, sig_s])

        return Sig_t, Sig_f, Sig_g, Sig_s
    


    
    def get_macro_cross_sections_groupwise(self):
        Sig_f = 1 
        Sig_g = 1.4
        Sig_s = 1 
        Sig_t = Sig_f + Sig_g + Sig_s
        return Sig_t, Sig_f, Sig_g, Sig_s


