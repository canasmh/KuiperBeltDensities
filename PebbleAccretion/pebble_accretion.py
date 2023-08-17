import numpy as np
import scipy.special as sp
from constants import AU_TO_CM, G, M_SUN
from gas_properties import column_density, gas_temp, kep_frequency, scale_height


def s_function(r_acc, h_p):
    y = (r_acc / 2 / h_p) ** 2
    I0 = sp.iv(0, y) 
    I1 = sp.iv(1, y) 
    S = np.exp(-y) * (I0 + I1)
    if np.isnan(S) == True or np.isinf(S) == True: 
            S = 2 / np.sqrt(2 * np.pi * y)
 
    return S


def focus_accretion(mass, density, rho_d, Omega, St, deltav, H_g, alpha=1e-4):
    radius = (3 * mass / (4 * np.pi * density)) ** (1 / 3) 
    vesc = np.sqrt(2 * G * mass / radius)
    tau_f = St / Omega
    Stk = deltav * tau_f / radius

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


if __name__ == "__main__":
    from constants import M_PLUTO
    from pebbles import Pebbles
    import matplotlib.pyplot as plt

    R = 20 * AU_TO_CM
    T = gas_temp(R)
    H = scale_height(T, R)
    rho_g = column_density(R) / (np.sqrt(2 * np.pi) * H)
    Omega = kep_frequency(R)
    t_orb = 2 * np.pi / Omega
    Z = 0.04
    n_pla = 100
    n_pebbles = 500
    planets = np.logspace(-5, 5, n_pla) * M_PLUTO
    pebbles = Pebbles(n_pebbles=n_pebbles, a_min=1e-4, a_max=1.0, rho_sil=3.5, model='bimodal', rho_ice=1.0, scale_height=H, gas_density=rho_g)
    fa_da = pebbles.volume_density_distribution(rho_d=rho_g * Z, alpha=1e-4)

    focus = np.zeros([n_pla, n_pebbles])
    bondi = np.zeros([n_pla, n_pebbles])
    hill = np.zeros([n_pla, n_pebbles])
    actual = np.zeros([n_pla, n_pebbles])

    transition_mass_da = np.sqrt(1. / 3) * 3000 ** 3 / G / Omega * (1 / 8) / pebbles.St

    for i in range(n_pla):
        for j in range(n_pebbles):
            focus[i, j] = focus_accretion(planets[i], 3.0, fa_da[j], Omega, pebbles.St[j], 3000, H, alpha=1e-4)
            bondi[i, j] = bondi_accretion(planets[i], R, fa_da[j], Omega, pebbles.St[j], 3000, H, chi=0.4, gamma=0.65, alpha=1e-4)
            hill[i, j] = hill_accretion(planets[i], R, fa_da[j], Omega, pebbles.St[j], 3000, H, chi=0.4, gamma=0.65, alpha=1e-4)

            if planets[i] <= transition_mass_da[j]:
                actual[i, j] = max(focus[i, j], bondi[i, j])

            else:
                actual[i, j] = hill[i, j]
    
    plt.plot(planets / M_PLUTO, np.sum(focus, axis=1) * t_orb / M_PLUTO, c='green', label="Focus Accretion")
    plt.plot(planets / M_PLUTO, np.sum(bondi, axis=1) * t_orb / M_PLUTO, c='magenta', label="Bondi Accretion")
    plt.plot(planets / M_PLUTO, np.sum(hill, axis=1) * t_orb / M_PLUTO, c='orange', label="Hill Accretion")
    plt.plot(planets / M_PLUTO, np.sum(actual, axis=1) * t_orb / M_PLUTO, c='k', label="Actual Accretion")
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'M (M$_{\rm{Pluto}}$)')
    plt.ylabel(r'$\dot{\rm{M}}$ (M$_{\rm{Pluto}}$/Orbit)')
    plt.xlim(1e-5, 1e5)
    plt.ylim(1e-12, 1e2)
    plt.legend()
    plt.show()