import numpy as np
from constants import *

def sound_speed(T):
    result = np.sqrt(T * C_P * (GAMMA - 1))

    return result


def kep_frequency(R):
    result = np.sqrt(G * M_SUN / R ** 3)

    return result


def scale_height(T, R):

    result = sound_speed(T) / kep_frequency(R)

    return result


def column_density(r, m_disk=0.04*M_SUN, r_c=50*AU_TO_CM):

    return m_disk / (2 * np.pi * r_c ** 2) * (r / r_c) ** (-1) * np.exp(- r / r_c)


def gas_temp(r):

    return 150 *  (r / AU_TO_CM) ** (-3 / 7)


def toomre_q(c_s, omega, sigma):

    return c_s * omega / (np.pi * G * sigma)


if __name__ == "__main__":
    R = [1, 10, 20, 30, 40, 50]
    r_c = 50 * AU_TO_CM

    for r in R:
        r *= AU_TO_CM
        
        T = gas_temp(r)
        print(f"R: {r / AU_TO_CM } AU")
        print(f"T: {T} K")
        print(f"Sound speed: {sound_speed(T): .3e} cm/s")
        print(f"Orbital Frequency: {kep_frequency(r): .3e} 1/s")
        print(f"Scale Height: {scale_height(T, r): .3e} cm")
        print(f"Column Density: {column_density(r): .3e} g/cm2")
        print(f"Column Density (non-truncated): {column_density(r) / np.exp(- r / r_c): .3e} g/cm2")
        print(f"Toomre Q (non-truncated vs. truncated): {toomre_q(sound_speed(T), kep_frequency(r), column_density(r) / np.exp(- r / r_c)): .2f} - {toomre_q(sound_speed(T), kep_frequency(r), column_density(r)): .2f}")
        print("")
