from pandas import read_csv as csv
import numpy as np


class StreamingInstabilityData:

    def __init__(self, rho_ice, rho_sil, unit_mass, initial_porosity=0.5, file_path="./data/si-data.csv", **kwargs):

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

        if "add_mass" in kwargs.keys():
            
            required_keys = ["n_bins", "m_per_bin", "min_mass", "max_mass", "min_dens", "max_dens"]
            for key in required_keys:
                if key in kwargs.keys():
                    continue
                else: 
                    self.add_mass = False
                    raise Exception(f"For additional masses, the following inputs need to be defined: {required_keys}")
            
            if kwargs["n_bins"] <= 0:
                raise ValueError("n_bins must be greater than 0")
            if kwargs["m_per_bin"] <= 0:
                raise ValueError("m_per_bin must be greater than 0")
            
            if kwargs["min_mass"] > kwargs["max_mass"]:
                raise ValueError("min_mass must be greater than max_mass")

            if kwargs["min_dens"] > kwargs["max_dens"]:
                raise ValueError("min_dens must be greater than max_dens")

            
            self.n_bins = int(kwargs["n_bins"])
            self.m_per_bin = int(kwargs["m_per_bin"])
            self.min_mass = float(kwargs["min_mass"])
            self.max_mass = float(kwargs["max_mass"])
            self.min_dens = float(kwargs["min_dens"])
            self.max_dens = float(kwargs["max_dens"])
            
            self.add_mass = True
        
        else:
            self.add_mass = False

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

        self.ice_fraction = [n_ice[i] / (n_ice[i] + n_sil[i]) for i in range(self.n_mass)]
        self.mass = [(n_ice[i] + n_sil[i]) * mpar for i in range(self.n_mass)]
        self.density = [self.porosity[i] * ((self.rho_ice * self.ice_fraction[i]) + (self.rho_sil * (1 - self.ice_fraction[i]))) for i in range(self.n_mass)]

        if self.add_mass:
            masses_to_add = np.linspace(self.min_mass, self.max_mass, self.n_bins)
            for _ in range(self.m_per_bin):
                for j in range(self.n_bins):
                    new_mass = masses_to_add[j]
                    new_density = self.min_dens + np.random.rand() * (self.max_dens - self.min_dens)
                    new_ice_fraction = (new_density - self.rho_sil * self.initial_porosity) / ((self.rho_ice - self.rho_sil) * self.initial_porosity)
                    self.mass += [new_mass]
                    self.density += [new_density]
                    self.ice_fraction += [new_ice_fraction]
        
        self.ice_fraction = np.array(self.ice_fraction)
        self.mass = np.array(self.mass)
        self.density = np.array(self.density)
        self.n_mass += (self.m_per_bin * self.n_bins)
    


if __name__ == "__main__":
    kbos = StreamingInstabilityData(rho_ice=1, rho_sil=3.5, unit_mass=2.823973078884959e+28, add_mass=True, min_dens=0.5, max_dens=0.6, min_mass=1e22, max_mass=1e23, n_bins=10, m_per_bin=30)
    # print(kbos.mass)
    # print(kbos.density)
    # print(kbos.ice_fraction)
    