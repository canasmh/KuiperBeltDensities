import numpy as np


class Pebbles:

    def __init__(self, n_pebbles, a_sil, a_ice, rho_sil, model, rho_ice=None):
        
        if model != "decreasing" or model != "constant" or model != "bimodal":
            raise ValueError("Model must be one of following: decreasing, constant, bimodal")

        elif rho_ice is None and model != "constant":
            raise ValueError("You must specify rho_ice if not using 'constant' model")

        else:
            self.a_sil = a_sil
            self.a_ice = a_ice
            self.rho_sil = rho_sil
            self.rho_ice = rho_ice
            self.model = model
            self.n_pebbles = n_pebbles
    
    def __get_pebbles_radius(self):
        self.raidus = np.logspace(np.log10(self.a_sil), np.log10(self.a_ice), self.n_pebbles)