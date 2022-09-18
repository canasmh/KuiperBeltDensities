from constants import AU_TO_CM, M_PLUTO
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

# Pebble Distribution
pebbles = Pebbles(n_pebbles=100, a_min=1e-4, a_max=1.0, rho_sil=3.5, model='bimodal', rho_ice=1.0, scale_height=H, gas_density=rho_g)
f_da = pebbles.volume_density_distribution(rho_d, alpha=1e-4)
W_da = pebbles.column_density_distribution(Z=Z, gas_column_density=sigma_g)

plt.plot(pebbles.radius, f_da)
plt.xscale('log')
plt.yscale('log')
plt.show()
    