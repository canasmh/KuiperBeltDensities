from pandas import read_csv as csv
import numpy as np


class StreamingInstabilityData:

    def __init__(self, rho_ice, rho_sil, initial_porosity=0.5, file_path="./data/si-data.csv"):
        self.rho_ice = rho_ice
        self.rho_sil = rho_sil
        self.file_path = file_path
        self.initial_porosity = initial_porosity
        self.__get_planetesimal_properties()

    def __read_data(self):
        data = csv(self.file_path)
        
        return data.n_ice, data.n_sil

    def __get_planetesimal_properties(self):
        n_ice, n_sil = self.__read_data()
        self.porosity = np.repeat(self.initial_porosity, len(n_ice))


if __name__ == "__main__":
    kbos = StreamingInstabilityData(rho_ice=1, rho_sil=3.5)
    
