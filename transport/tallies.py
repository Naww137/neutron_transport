
import numpy as np



class tallies:

    def __init__(self, 
                Emin=None, Emax=None, iEbins=None,
                Xmin=None, Xmax=None, iXbins=None):

        ### setup energy bins
        if iEbins is not None:
            self.setup_Ebins(Emin, Emax, iEbins)
            # self.collision_rate_tally = []

        ### setup spatial bins
        if iXbins is not None:
            self.setup_Xbins(Xmin,Xmax, iXbins)
            
        ### init runtime variables
        self.simulated_neutrons = 0
        self.generations = 0
        self.k_estimates = []
        self.population_variance = []
        self.estimator_variance = []

        return

    def setup_Ebins(self, Emin, Emax, iEbins):
        vEbins = np.logspace(np.log10(Emin),np.log10(Emax),iEbins+1) #edges
        dEbins = np.zeros(len(vEbins)-1)
        Ebins = np.zeros(len(vEbins)-1) # midpoints
        for i in range(1,len(vEbins)):
            dEbins[i-1] = vEbins[i] - vEbins[i-1] 
            Ebins[i-1] = (vEbins[i]+vEbins[i-1])/2 
        self.Emin = Emin
        self.Emax = Emax
        self.iEbins = iEbins
        self.vEbins = vEbins
        self.dEbins = dEbins
        self.Ebins = Ebins # midpoints
        self.Etally = np.zeros((2,len(Ebins))) # tally container
        return
    
    def setup_Xbins(self, Xmin, Xmax, iXbins):
        vXbins = np.linspace(Xmin,Xmax,iXbins+1) #edges
        dXbins = np.zeros(len(vXbins)-1)
        Xbins = np.zeros(len(vXbins)-1) # midpoints
        for i in range(1,len(vXbins)):
            dXbins[i-1] = vXbins[i] - vXbins[i-1] 
            Xbins[i-1] = (vXbins[i]+vXbins[i-1])/2 
        self.Xmin = Xmin
        self.Xmax = Xmax
        self.iXbins = iXbins
        self.vXbins = vXbins
        self.dXbins = dXbins
        self.Xbins = Xbins # midpoints
        self.Xtally = np.zeros((2,len(Xbins))) # tally container

        self.MCsource = np.zeros((2,len(Xbins)))
        self.FMsource = np.zeros((2,len(Xbins)))
        return
    

    def tally_source(self, source_particle_location):
        bindex = np.searchsorted(self.vXbins, source_particle_location)-1
        
        return

    # def setup_ktally(self, Emin, Emax):
    #     self.intragen_ktally = np.zeros((2,1))
    #     return
    
    def tally_energy(self, E, sigma):
        bindex = np.searchsorted(self.vEbins, E)-1
        self.Etally[0, bindex] += 1/sigma/self.dEbins[bindex]
        self.Etally[1, bindex] += (1/sigma/self.dEbins[bindex])**2
        return
        
    def save_generation_tally(self, N):
        ### calculate mean k estimate, population variance estimate, and mean k estimate variance
        k = self.first_moment/N
        population_variance = (N/(N-1)) * (((self.second_moment)/N) - ((self.first_moment/N)**2)) 
        estimator_variance = (1/np.sqrt(N)) * np.sqrt(self.population_variance)
        ### save stuff
        self.k_estimates.append(k)
        self.population_variance.append(population_variance)
        self.estimator_variance.append(estimator_variance)
        self.simulated_neutrons += N
        self.generations += 1
        return
    
    def calculate_flux(self, Etally, N):
        mean = Etally[0]/N
        variance = 1/np.sqrt(N) * (Etally[1]/N - mean**2)
        collision_based = np.array([mean, variance])
        return collision_based
    
    def final_analysis(self):
        final_k_estimate = np.mean(self.k_estimates)
        final_population_variance = (1/(self.generations-1)) * np.sum((self.k_estimates-final_k_estimate)**2)
        final_estimator_variance =  (1/np.sqrt(self.generations))*np.sqrt(final_population_variance) 
        collision_based_scalar_flux = self.calculate_flux(self.Etally, self.simulated_neutrons)
        return final_k_estimate, final_estimator_variance, collision_based_scalar_flux


    def reset_generation_tally(self):
        self.first_moment = 0 # neutrons
        self.second_moment = 0 # neutrons^2
        return


# class intergen(tallies):

#     def __init__(self, Emin,Emax,iEbins):
#         tallies.__init__(self, Emin,Emax,iEbins)
#         self.simulated_neutrons = 0
#         self.generations = 0
#         self.k_estimates = []
#         self.population_variance = []
#         self.estimator_variance = []

#     def tally_generation(self, gtal, N):

#         self.k_estimates.append(gtal.k)
#         self.population_variance.append(gtal.population_variance)
#         self.estimator_variance.append(gtal.estimator_variance)

#         self.simulated_neutrons += N
#         self.generations += 1

#         return
    
#     def final_analysis(self):
#         final_k_estimate = np.mean(self.k_estimates)
#         final_population_variance = (1/(self.generations-1)) * np.sum((self.k_estimates-final_k_estimate)**2)
#         final_estimator_variance =  (1/np.sqrt(self.generations))*np.sqrt(final_population_variance) 
#         return final_k_estimate, final_estimator_variance


