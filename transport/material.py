
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

        return sig_t, sig_g, sig_s, sig_f

    def get_macro_cross_sections(self, energy):

        sig_t, sig_g, sig_s, sig_f = self.get_micro_cross_sections(energy)
        Sig_t, Sig_g, Sig_s, Sig_f = self.number_density*np.array([sig_t, sig_g, sig_s, sig_f])

        return Sig_t, Sig_g, Sig_s, Sig_f
    


class material_MG:

    def __init__(self, Sig_f, Sig_g, Sig_s):

        self.Sig_f = Sig_f
        self.Sig_g = Sig_g
        self.Sig_s = Sig_s
        self.Sig_t = self.Sig_f + self.Sig_g + self.Sig_s

        return
    
    def get_macro_cross_sections_groupwise(self):

        return self.Sig_t, self.Sig_g, self.Sig_s, self.Sig_f


