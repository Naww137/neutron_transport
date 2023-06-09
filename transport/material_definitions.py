
from transport import functions_for_scattering_theory as xs
import numpy as np

class isotope:

    def __init__(self, pair_constants, resonance_ladder, number_density, nubar):

        self.pair_constants = pair_constants
        self.resonance_ladder = resonance_ladder
        self.number_density = number_density
        self.nubar = nubar

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
    

    def get_macro_cross_sections_groupwise(self):
        Sig_f = 1 
        Sig_g = 1.4
        Sig_s = 1 
        Sig_t = Sig_f + Sig_g + Sig_s
        return Sig_t, Sig_g, Sig_s, Sig_f


class material:

    # perhaps add region definition
    def __init__(self, isotopes, nubar = 2.88,
                                constant_scattering=0):
        self.isotopes = isotopes
        self.constant_scattering = constant_scattering
        self.nubar = nubar
        return
    
    def get_macro_cross_sections(self, energy):
        m_Sig_f, m_Sig_g, m_Sig_s = 0, 0, 0
        for isotope in self.isotopes:
            _, i_sig_g, i_sig_s, i_sig_f = isotope.get_micro_cross_sections(energy)
            i_Sig_g, i_Sig_s, i_Sig_f = isotope.number_density*np.array([i_sig_g, i_sig_s, i_sig_f])
            m_Sig_f += i_Sig_f
            m_Sig_g += i_Sig_g
            m_Sig_s += i_Sig_s
        m_Sig_s += self.constant_scattering
        m_Sig_t  = m_Sig_f + m_Sig_g + m_Sig_s
        return m_Sig_t, m_Sig_g, m_Sig_s, m_Sig_f


class material_MG:

    def __init__(self, nubar, Sig_f, Sig_g, Sig_s):

        self.nubar = nubar
        self.Sig_f = Sig_f
        self.Sig_g = Sig_g
        self.Sig_s = Sig_s
        self.Sig_t = self.Sig_f + self.Sig_g + self.Sig_s

        return
    
    def get_macro_cross_sections(self):

        return self.Sig_t, self.Sig_g, self.Sig_s, self.Sig_f
