import numpy as np


class Pebbles:
    # TODO: Need to get ice fraction pebbles. also need to formulate the volume and column density distribution

    def __init__(self, n_pebbles, a_min, a_max, rho_sil, model, rho_ice=None):
        
        if model != "decreasing" and model != "constant" and model != "bimodal":
            raise ValueError("Model must be one of following: decreasing, constant, bimodal")

        elif rho_ice is None and model != "constant":
            raise ValueError("You must specify rho_ice if not using 'constant' model")

        if model != "constant":
            if rho_ice >= rho_sil:
                raise ValueError("Ice density must be larger than silicate density")

        if a_min > a_max:
            raise ValueError("a_min must be greater than a_max")

        self.a_sil = a_min
        self.a_ice = a_max
        self.rho_sil = rho_sil
        self.rho_ice = rho_ice
        self.model = model
        self.n_pebbles = n_pebbles
        self.bimodal_split = 0.1  # Bimodal distribution splits at 1 mm
        self.__get_pebble_radius()
        if model == "constant":
            self.density = np.repeat(rho_sil, n_pebbles)
        else:
            self.__get_pebble_density()

    def __get_pebble_radius(self):
        self.radius = np.logspace(np.log10(self.a_sil), np.log10(self.a_ice), self.n_pebbles)

    def __get_pebble_density(self):

        rhops_test = 0
        q = 0
        
        if self.model == "decreasing":
            while round(rhops_test, 3) != self.rho_sil and q < 1:
                q += 1e-5
                rhops = self.rho_ice * (self.radius / self.radius.max()) ** (-q)
                rhops_test = rhops.max()

        else:
            rhops = np.repeat(0.0, self.n_pebbles)
            for i in range(self.n_pebbles):
                if self.radius[i] <= self.bimodal_split:
                    continue
                else:
                    break
            
            while round(rhops_test, 3) != self.rho_sil and q < 1:
                q += 1e-5
                rhops[i:] = self.rho_ice * (self.radius[i:] / self.radius.max()) ** (-q)
                rhops_test = rhops[i:].max()

            rhops[:i] = np.repeat(self.rho_sil, i)
        
        self.rhops = rhops
        self.q = q


if __name__ == "__main__":
    pebbles = Pebbles(n_pebbles=100, a_min=1e-4, a_max=1, rho_sil=3.7, rho_ice=0.5, model="bimodal")