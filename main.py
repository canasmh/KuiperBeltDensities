from accretion_integrator import rk4
from constants import AU_TO_CM, M_PLUTO, M_SUN, YRS_TO_SEC, G
from gas_properties import sound_speed, scale_height, kep_frequency, gas_temp, column_density, toomre_q
from pebbles import Pebbles
from plotter import plot_kbo_data, animate_figures
from reader import StreamingInstabilityData
import numpy as np
import matplotlib.pyplot as plt

R = 20 * AU_TO_CM
NOTES = "m_truncated_CO"
# Disk Conditions
T = max(gas_temp(R), 30)
r_c = 50 * AU_TO_CM
sigma_g =  column_density(R, m_disk=0.1 * M_SUN)  # FOR CO ICE USE M_DISK = 0.05 * M_SUN # For constant use 0.1
c_s = sound_speed(T)
H = scale_height(T, R)
Omega = kep_frequency(R)
Q = toomre_q(c_s, Omega, sigma_g)

if Q < 10:
    print(f"WARNING: Disk is not non self-gravitating -- Q = {Q}")
    if round(Q) <= 1:
        raise ValueError(f"Disk is unstable: Toomre Q is {Q:.3f}")
else:
    print(f"Toomre Q: {Q}")

rho_g_init = sigma_g / (np.sqrt(2 * np.pi) * H)
rho_g = rho_g_init
Z = 0.01 # Dust to gas ratio
rho_d = Z * rho_g

# Streaming Instability results
rho_sil = 3.0
kbos = StreamingInstabilityData(rho_ice=0.89, rho_sil=rho_sil, unit_mass=rho_g * (H ** 3))
kbos.add_masses(n_bins=50, m_per_bin=4, min_dens=min(kbos.density), max_dens=min(kbos.density) + 0.1, min_mass=5e-4 * M_PLUTO, max_mass=1e-2 * M_PLUTO)
n_mass = kbos.n_mass

# Pebble Distribution
n_pebbles = 50
model='constant'
pebbles = Pebbles(n_pebbles=n_pebbles, a_min=1e-4, a_max=1.0, rho_ice=0.89, rho_sil=rho_sil, model=model, scale_height=H, gas_density=rho_g)
f_da = pebbles.volume_density_distribution(rho_d, alpha=1e-4)
W_da = pebbles.column_density_distribution(Z=Z, gas_column_density=sigma_g)
delta_v = 3000  # sub-keplerian velocity
transition_mass = np.sqrt(1 / 3) * delta_v ** 3 / G / Omega * (1 / 8) / pebbles.St # If mass > than this value, its in the Hill regime

# Run parameters
IVAR = 12  # The last printed 'var' file from our simulation. Var prints every half orbit
t_orb = 2 * np.pi / Omega
tol = 0.01  # The maximum mass accretion rate per dt will be 0.5% of the mass of corresponding KBO
dt_init = 100 * t_orb  # Set the initial time step to be 100 orbits
dt = dt_init
t_min = IVAR * np.pi / (2 * np.pi) * t_orb  # Starting time is 6 orbits
t_max = 5.8271e6 * YRS_TO_SEC # For CO ice use 5.728e6 * YRS_TO_SEC  # Maximum time of simulation
# 5.81e6 is to low for silicates..
# 5.83e6 is a little too high for silicates...
t = t_min  # Set time to t_min
it = 0
icadence = 3
m_max = 1e1 * M_PLUTO

# Print output headers
print(f"{'it':=^10}{'t(yrs)':=^15}{'dt (yrs)':=^15}{'max(m) Pluto':=^20}{'rho(max_m)':=^15}{'ice % (max_m)':=^15}{'phi(max_m)':=^15}{'dm(max_m)':=^15}{'ice % added (max_m)':=^15}{'optimal a (max_m)':=^20}")
while t < t_max and max(kbos.mass) < m_max:
    
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
        radius_to_fully_compact = 1500 * 1e5

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

    if it % icadence == 0:
        imax = np.where(kbos.mass == max(kbos.mass))
        dm_max = float(np.sum(M_ADDED[:,imax]))
        idm_max = np.where(M_ADDED[:,imax] == np.max(M_ADDED[:, imax]))[0][0]
        plot_kbo_data(t=t / YRS_TO_SEC, i=it // icadence, si_data=kbos, savefig=True)
        # print(f"{it:^10}{t / YRS_TO_SEC:^15.3e}{dt / YRS_TO_SEC:^15.3e}{kbos.mass[imax][0] / M_PLUTO:^20.5e}{kbos.density[imax][0]:^15.5f}{kbos.ice_fraction[imax][0] * 100:^15.3f}{kbos.porosity[imax][0]:^15.3f}{dm_max / (dt / t_orb) / M_PLUTO:^15.3e}{ice_frac_added_list[imax[0][0]] * 100:^15.3f}")#{pebbles.radius[idm_max]:^20.3e}")

        print(f"{it:^10}{t / YRS_TO_SEC:^15.3e}{dt / YRS_TO_SEC:^15.3e}{kbos.mass[imax][0] / M_PLUTO:^20.5e}{kbos.density[imax][0]:^15.5f}{kbos.ice_fraction[imax][0] * 100:^15.3f}{kbos.porosity[imax][0]:^15.3f}{dm_max / (dt / t_orb) / M_PLUTO:^15.3e}{ice_frac_added_list[int(imax[0][0])] * 100:^15.3f}{pebbles.radius[idm_max]:^20.3e}")
        

    t += dt
    it += 1

if model == "bimodal":
    model_alias = "bi"

elif model == "constant":
    model_alias = "k"
    
else:
    model_alias = "dec"


filename = f"R_{int(R / AU_TO_CM)}_{model_alias}_rhopsmax_{'_'.join(str(rho_sil).split('.'))}_tmax_{int(t_max / YRS_TO_SEC)}_{NOTES}"
animate_figures(filename=filename)

np.savez(filename, mass=kbos.mass, density=kbos.density, ice_fraction=kbos.ice_fraction, porosity=kbos.porosity)
    