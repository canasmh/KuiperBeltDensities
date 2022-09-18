import numpy as np
import scipy.special as sp
from constants import G, M_SUN


def s_function(r_acc, h_p):
    y = (r_acc / 2 / h_p) ** 2
    I0 = sp.iv(0, y) 
    I1 = sp.iv(1, y) 
    S = np.exp(-y) * (I0 + I1)
    if np.isnan(S) == True or np.isinf(S) == True: 
            S = 2 / np.sqrt(2 * np.pi * y)
 
    return S


def focus_accretion(mass, density, rho_d, Omega, St, deltav, H_g, alpha=1e-4):
    tau_f = St / Omega
    Stk = deltav * tau_f / radius
    radius = (3 * mass / (4 * np.pi * density)) ** (1 / 3) 
    vesc = np.sqrt(2 * G * mass / radius)

    if (Stk <= 1):
        return 0.
        
    else:
        H_p = H_g / np.sqrt(1 + St / alpha)
        
        S = s_function(radius, H_p)
         
        return np.pi * radius ** 2 * rho_d * S * deltav * (1 + vesc ** 2 / deltav ** 2) 


def bondi_accretion(m, r, rho_d, Omega, St, deltav, H_g, chi=0.4, gamma=0.65, alpha=1e-4):

    # Bondi radius and Bondi time scale.
    r_bondi = G * m / deltav ** 2
    t_bondi = r_bondi / deltav
    R_Hill = r * (1 / 3 * m / M_SUN) ** (1 / 3)
    tp = G * m / (deltav + Omega * R_Hill) ** 3
    t_f = St / Omega
    tf_over_tbondi = t_f / t_bondi
    r_acc_hat = 2 * np.sqrt(tf_over_tbondi) * r_bondi
    fac = np.exp(-chi * (t_f / tp) ** gamma)
    Racc_Bondi = r_acc_hat * fac
    H_p = H_g / np.sqrt(1 + St / alpha)
    S = s_function(Racc_Bondi, H_p)
    dv = deltav + Omega * Racc_Bondi

    return np.pi * Racc_Bondi ** 2 * rho_d * dv * S

def hill_accretion(m, r, rho_d, Omega, St, deltav, H_g, chi=0.4, gamma=0.65, alpha=1e-4):

    chi = 0.4
    gamma = 0.65
    
    tau_f = St / Omega
    H_p = H_g / np.sqrt(1 + St / alpha)
    r_Hill = r * (1 / 3 * m / M_SUN) ** (1 / 3)
    tp = G * m / (deltav + Omega * r_Hill) ** 3
    fac = np.exp(-chi * (tau_f / tp) ** gamma)
    
    Racc_Hill = np.cbrt(St/ 0.1) * r_Hill * fac
    S = s_function(Racc_Hill, H_p)

    dv_hill = deltav + Omega * Racc_Hill

    return np.pi *  Racc_Hill ** 2 * S * rho_d * dv_hill

