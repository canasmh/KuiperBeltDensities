from accretion_integrator import rk4
from constants import AU_TO_CM, M_PLUTO, YRS_TO_SEC, G
from gas_properties import sound_speed, scale_height, kep_frequency, gas_temp, column_density
from pebbles import Pebbles
from reader import StreamingInstabilityData
import numpy as np
import matplotlib.pyplot as plt

R = 20 * AU_TO_CM

# Disk Conditions
T = max(gas_temp(R), 30)
sigma_g =  column_density(R)
c_s = sound_speed(T)
H = scale_height(T, R)
Omega = kep_frequency(R)
rho_g_init = sigma_g / (np.sqrt(2 * np.pi) * H)
rho_g = rho_g_init
Z = 0.04 # Dust to gas ratio
rho_d = Z * rho_g

# Streaming Instability results
kbos = StreamingInstabilityData(rho_ice=1, rho_sil=3.0, unit_mass=rho_g * (H ** 3))
kbos.add_masses(n_bins=50, m_per_bin=3, min_dens=min(kbos.density), max_dens=min(kbos.density) + 0.1, min_mass=1e-3 * M_PLUTO, max_mass=1e-2 * M_PLUTO)
n_mass = kbos.n_mass

# Pebble Distribution
n_pebbles = 100
pebbles = Pebbles(n_pebbles=n_pebbles, a_min=1e-4, a_max=1.0, rho_sil=3.5, model='bimodal', rho_ice=1.0, scale_height=H, gas_density=rho_g)
f_da = pebbles.volume_density_distribution(rho_d, alpha=1e-4)
W_da = pebbles.column_density_distribution(Z=Z, gas_column_density=sigma_g)
delta_v = 3000  # sub-keplerian velocity
transition_mass = np.sqrt(1 / 3) * delta_v ** 3 / G / Omega * (1 / 8) / pebbles.St # If mass > than this value, its in the Hill regime

# Run parameters
IVAR = 12  # The last printed 'var' file from our simulation. Var prints every half orbit
t_orb = 2 * np.pi / Omega
tol = 0.005  # The maximum mass accretion rate per dt will be 0.5% of the mass of corresponding KBO
dt_init = 10 * t_orb  # Set the initial time step to be 10 orbits
t_min = IVAR * np.pi / (2 * np.pi) * t_orb  # Starting time is 6 orbits
t_max = 4.5e6 * YRS_TO_SEC  # Simulation will run for 4.5 Million years
t = t_min  # Set time to t_min
it = 0

while t < t_max:
    dt = dt_init
    rho_g = rho_g_init * np.exp(-t / t_max)
    pebbles.stokes_number(H, rho_g)
    transition_mass = np.sqrt(1 / 3) * delta_v ** 3 / G / Omega * (1 / 8) / pebbles.St
    f_da = pebbles.volume_density_distribution(rho_d, alpha=1e-4)

    M_BONDI = np.zeros([n_pebbles, n_mass])
    M_FOCUS = np.zeros([n_pebbles, n_mass])
    M_ADDED = np.zeros([n_pebbles, n_mass])

    # Get the time step
    mass_exceeded = True

    while mass_exceeded:
        for m in range(n_mass):
            for p in range(n_pebbles):
                if kbos.mass[m] < transition_mass[p]:
                    M_FOCUS[p, m] = rk4(dt, acc_type='focus', m=kbos.mass[m], rho=kbos.density[m], rho_d=f_da[p], Omega=Omega, St=pebbles.St[p], deltav=delta_v, H_g=H)
                    M_BONDI[p, m] = rk4(dt, acc_type='bondi', m=kbos.mass[m], r=R, rho_d=f_da[p], Omega=Omega, St=pebbles.St[p], deltav=delta_v, H_g=H)
                    M_ADDED[p, m] = max(M_FOCUS[p, m], M_BONDI[p, m])

                else:
                    M_ADDED[p, m] = rk4(dt, acc_type='hill', m=kbos.mass[m], r=R, rho_d=f_da[p], Omega=Omega, St=pebbles.St[p], deltav=delta_v, H_g=H)
            
            if np.sum(M_ADDED[:, m]) > tol * kbos.mass[m]:
                dt *= 0.9
                break

            elif m == n_mass - 1:
                mass_exceeded = False

    ice_frac_added_list = []
    for m in range (n_mass):
        density_added, ice_frac_added, ice_frac_for_list = (0, 0, 0)
        old_mass = kbos.mass[m]
        kbos.mass[m] += np.sum(M_ADDED[:, m])

        for p in range(n_pebbles):
            dm = M_ADDED[p, m]
            f_mass_rhops = dm / kbos.mass[m]
            density_added += (pebbles.density[p] * f_mass_rhops)
            ice_frac_added += (pebbles.ice_fraction[p] * f_mass_rhops)
            ice_frac_for_list += pebbles.ice_fraction[p] * dm / np.sum(M_ADDED[:, m])
        ice_frac_added_list.append(ice_frac_for_list)

        kbos.density[m] = kbos.density[m] * (old_mass / kbos.mass[m]) + density_added * (1 - kbos.porosity[m])
        kbos.ice_fraction[m] = kbos.ice_fraction[m] * (old_mass / kbos.mass[m]) + ice_frac_added

        # Radius a KBO must have to be fully compact
        radius_to_fully_compact = 1455 * 1e5

        # The porosity of kbo, based on its size
        new_porosity = min(np.log10(radius_to_fully_compact) - np.log10(kbos.radius(i=m)), 0.5)

        if new_porosity < 0:
            new_porosity = 0

        # Remove the porosity (i.e., make it fully compact)
        kbos.density[m] /= (1 - kbos.porosity[m])

        # Replace porosity with new porosity
        kbos.porosity[m] = new_porosity

        # Get the new density
        kbos.density[m] *= (1 - kbos.porosity[m])

    t += dt
    