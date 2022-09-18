import numpy as np


class Pebbles:
    # TODO: Need to formulate the volume and column density distribution

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
            self.q = 0
            self.density = np.repeat(rho_sil, n_pebbles)
            self.ice_fraction = np.repeat(0.0, n_pebbles)

        else:
            self.__get_pebble_density()
        
        self.__get_pebble_mass()

    def __get_pebble_radius(self):
        self.radius = np.logspace(np.log10(self.a_sil), np.log10(self.a_ice), self.n_pebbles)

    def __get_pebble_density(self):

        rhops_test = 0
        q = 0
        
        if self.model == "decreasing":
            while round(rhops_test, 3) != self.rho_sil:
                q += 1e-5
                rhops = self.rho_ice * (self.radius / self.radius.max()) ** (-q)
                rhops_test = rhops.max()

                if rhops_test > self.rho_sil:
                    raise ValueError("Could not get the correct power law. Please ensure you  are using less than 4 significant figures.")

        else:
            rhops = np.repeat(0.0, self.n_pebbles)
            for i in range(self.n_pebbles):
                if self.radius[i] <= self.bimodal_split:
                    continue
                else:
                    break
            
            while round(rhops_test, 3) != self.rho_sil:
                q += 1e-5
                rhops[i:] = self.rho_ice * (self.radius[i:] / self.radius.max()) ** (-q)
                rhops_test = rhops[i:].max()

                if rhops_test > self.rho_sil:
                    raise ValueError("Could not get the correct power law. Please ensure you  are using less than 4 significant figures.")

            rhops[:i] = np.repeat(self.rho_sil, i)
        
        self.density = rhops
        self.ice_fraction = abs(1 - ((self.density - self.density[-1]) / (self.density[0] - self.density[-1])))
        self.q = q
    
    def __get_pebble_mass(self):
        self.mass = self.density * 4 / 3 * np.pi * self.radius ** 3


if __name__ == "__main__":
    import time
    n_par = 10
    models = ["decreasing", "bimodal", "constant"]
    sil_densities = [2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 3.7896, 3.956]
    ice_densities = [0.5, 1.0]
    a_min = 1e-4
    a_max = 1

    for model in models:
        print(f"{model:-^120}")
        time.sleep(1)
        for sil_density in sil_densities:
            for ice_density in ice_densities:
                try:
                    pebbles = Pebbles(n_pebbles=n_par, a_min=a_min, a_max=a_max, rho_sil=sil_density, model=model, rho_ice=ice_density)
                except ValueError:
                    print(f"Could not converge for density:\nrho_ice = {ice_density}\nrho_sil = {sil_density}")
                time.sleep(1)
                print(f"{' radius (cm) ':=^30}{' mass (g) ':=^30}{' density (g/cm3) ':=^30}{' ice fraction (%) ':=^30}")
                for i in range(pebbles.n_pebbles):
                    print(f"{pebbles.radius[i]:^30.3e}{pebbles.mass[i]:^30.3e}{pebbles.density[i]:^30.3f}{pebbles.ice_fraction[i] * 100:^30.3f}")
                print("")

