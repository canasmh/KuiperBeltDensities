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
        self.n_mass = len(n_ice)
        self.porosity = np.repeat(self.initial_porosity, self.n_mass)
        mgas_code = self.total_density * self.dx * self.dy * self.dz  # Sum of the gas mass in our simulation (code units)
        mpar_code = mgas_code * self.eps / self.npar
        mpar = mpar_code * self.unit_mass

        self.ice_fraction = np.array([n_ice[i] / (n_ice[i] + n_sil[i]) for i in range(self.n_mass)])
        self.mass = np.array([(self.n_ice[i] + self.n_sil[i]) * mpar for i in range(self.n_mass)])
        self.density = np.array([self.porosity[i] * ((self.rho_ice * self.ice_fraction[i]) + (self.rho_sil * (1 - self.ice_fraction[i]))) for i in range(self.n_mass)])


if __name__ == "__main__":
    kbos = StreamingInstabilityData(rho_ice=1, rho_sil=3.5)
    
