import numpy as np
import scipy.special as sp


def s_function(r_acc, h_p):
    y = (r_acc / 2 / h_p) ** 2
    I0 = sp.iv(0, y) 
    I1 = sp.iv(1, y) 
    S = np.exp(-y) * (I0 + I1)
    if np.isnan(S) == True or np.isinf(S) == True: 
            S = 2 / np.sqrt(2 * np.pi * y)
 
    return S


def focus_accretion(mass, density, rho_d, Omega, St, vesc, deltav, H_g, alpha=1e-4):
    tau_f = St / Omega
    Stk = deltav * tau_f / radius
    radius = (3 * mass / (4 * np.pi * density)) ** (1 / 3) 

    if (Stk <= 1):
        return 0.
        
    else:
        H_p = H_g / np.sqrt(1 + St / alpha)
        
        S = s_function(radius, H_p)
         
        return np.pi * radius ** 2 * rho_d * S * deltav * (1 + vesc ** 2 / deltav ** 2) 

