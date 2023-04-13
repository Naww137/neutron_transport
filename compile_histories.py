#%%
import numpy as np

#%%
# flux = []
# keff = []

# for i in range(10000):
#     file = f"./transport_result_{i}.csv"
#     # file = f"/Users/noahwalton/Library/Mobile Documents/com~apple~CloudDocs/References/NE 697 Advanced Transport Methods/URR Project/histories/transport_result_{i}.csv"
#     try:
#         flux_dat = np.genfromtxt(file, delimiter=',', skip_header=1)
#         if len(flux_dat) > 400:
#             flux_dat = flux_dat[0:-1]
#         k_unc_dat = np.loadtxt(file, delimiter=',', max_rows=1, usecols=(0,1))
#         flux.append(flux_dat)
#         keff.append(k_unc_dat)
#     except:
#         pass


# flux = np.array(flux)
# keff = np.array(keff)

# np.save('./fluxes.npy', flux)
# np.save('./keffs.npy', keff)

#%%

fluxes = np.load("/Users/noahwalton/Library/Mobile Documents/com~apple~CloudDocs/References/NE 697 Advanced Transport Methods/URR Project/fluxes.npy")
keffs = np.load("/Users/noahwalton/Library/Mobile Documents/com~apple~CloudDocs/References/NE 697 Advanced Transport Methods/URR Project/keffs.npy")
# %%
