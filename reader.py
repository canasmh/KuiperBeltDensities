from pandas import read_csv as csv
import numpy as np


class StreamingInstabilityData:

    def __init__(self, rho_ice, rho_sil, unit_mass, initial_porosity=0.5, file_path="./data/si-data.csv"):
        self.rho_ice = rho_ice  # Density of icy pebbles
        self.rho_sil = rho_sil  # Density of silicate pebbles
        self.file_path = file_path  # File path to the CSV file
        self.unit_mass = unit_mass  # Actual gas mass contained in one cubic scale height
        self.initial_porosity = initial_porosity  # Initial porosity of KBOs

        # The following parameters are from the simulation -- required to convert to physical units
        self.npar = 1.536e7  # Number of particles used in the simulation
        self.eps = 0.03  # Dust to gas ratio
        self.total_density = 16749076.820152447  # Sum of the gas density in our simulation (code units)
        self.dx, self.dy, self.dz = (0.00078125, 0.00078125, 0.00078125)  # Grid spacing in out simulation
        self.__get_planetesimal_properties()

    def __read_data(self):
        data = csv(self.file_path)
        
        return data.n_ice, data.n_sil

    def __get_planetesimal_properties(self):
        n_ice, n_sil = self.__read_data()
        self.porosity = np.repeat(self.initial_porosity, len(n_ice))
        mgas_code = self.total_density * self.dx * self.dy * self.dz  # Sum of the gas mass in our simulation (code units)
        mgas = mgas_code * self.unit_mass  # Same as above but physical units
        mpar_code = mgas_code * self.eps / self.npar
        mpar = mpar_code * self.unit_mass

        kbo_mass = []
        kbo_density = []
        kbo_ice_fraction = []

        for i in range(N_MASS):
            kbo_ice_fraction.append(sinks_n_ice[i] / (sinks_n_ice[i] + sinks_n_sil[i]))
            kbo_mass.append((sinks_n_ice[i] + sinks_n_sil[i]) * mpar)
            kbo_density.append(POROSITY * ((dens_ice * kbo_ice_fraction[i]) + (dens_sil * (1 - kbo_ice_fraction[i]))))


if __name__ == "__main__":
    kbos = StreamingInstabilityData(rho_ice=1, rho_sil=3.5)
    
