import numpy as np

# BOTH MODELS TAKEN FROM Muraoka, R. [https://doi.org/10.1063/5.0122861]
# Models for mach disk location
def md_loc(eta_0):
    return 0.67*np.sqrt(eta_0)

# Models for mach disk diameter
def md_dia(eta_e):
    return 5/2*np.log10(eta_e) - 3/4