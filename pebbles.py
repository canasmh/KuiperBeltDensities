import numpy as np
from constants import AU_TO_CM


class Pebbles:
    # TODO: Need to formulate the volume and column density distribution

    def __init__(self, n_pebbles, a_min, a_max, rho_sil, model, scale_height, gas_density, rho_ice=None):
        
        if model != "decreasing" and model != "constant" and model != "bimodal":
            raise ValueError("Model must be one of following: decreasing, constant, bimodal")

        elif rho_ice is None and model != "constant":
            raise ValueError("You must specify rho_ice if not using 'constant' model")

        if model != "constant":
            if rho_ice >= rho_sil:
                raise ValueError("Ice density must be larger than silicate density")

        if a_min > a_max:
            raise ValueError("a_min must be less than a_max")

        self.a_sil = a_min
        self.a_ice = a_max
        self.rho_sil = rho_sil
        self.rho_ice = rho_ice
        self.model = model
        self.n_pebbles = n_pebbles
        self.bimodal_split = 0.01  # Bimodal distribution splits at 0.1 mm
        self.__get_pebble_radius()
        
        if model == "constant":
            self.q = 0
            self.density = np.repeat(rho_sil, n_pebbles)
            self.ice_fraction = np.repeat(0.0, n_pebbles)

        else:
            self.__get_pebble_density()
        
        self.__get_pebble_mass()
        self.stokes_number(scale_height, gas_density)

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
                    raise ValueError("Could not get the correct power law. Please ensure you are using less than 4 significant figures.")

            rhops[:i] = np.repeat(self.rho_sil, i)

        if q >= 0.5:
            raise ValueError("Negative densities encounteres (q >= 0.5. This is resulting from needing a large exponent to get from rho_sil to rho_ice")
        
        self.density = rhops
        self.ice_fraction = abs(1 - ((self.density - self.density[-1]) / (self.density[0] - self.density[-1])))
        self.q = q
    
    def __get_pebble_mass(self):
        self.mass = self.density * 4 / 3 * np.pi * self.radius ** 3

    def stokes_number(self, scale_height, gas_density):
        self.St = self.radius * self.density / (scale_height * gas_density)
        
    def column_density_distribution(self, Z, gas_column_density):
        p = 0.5 + self.q
        W = 3 * (1 - p) * Z * gas_column_density / (4 * np.pi * self.density[-1] * np.sqrt(self.radius[-1])) * self.radius ** (-3.5)

        return self.mass * W * np.gradient(self.radius)

    def volume_density_distribution(self, rho_d, alpha):
        p = 0.5 + self.q
        f = 3 * (1 - p) * rho_d / (4 * np.pi * self.density[-1]* np.sqrt(self.radius[-1])) * np.sqrt(1 + self.St / alpha) * self.radius ** (-3.5)

        return self.mass * f * np.gradient(self.radius)


if __name__ == "__main__":
    import time
    from gas_properties import gas_temp, column_density, scale_height
    R = 20 * AU_TO_CM

    T = gas_temp(R)
    sigma_g = column_density(R)
    H = scale_height(T, R)
    rho_g = sigma_g / (np.sqrt(2 * np.pi) * H)
    Z = 0.04
    rho_d = Z * rho_g
    n_par = 20
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
                    pebbles = Pebbles(n_pebbles=n_par, a_min=a_min, a_max=a_max, rho_sil=sil_density, model=model, rho_ice=ice_density, scale_height=H, gas_density=rho_g)
                    pebbles.stokes_number(H, rho_g)
                    f_da = pebbles.volume_density_distribution(rho_d=rho_d, alpha=1e-4)
                    W_da = pebbles.column_density_distribution(Z, sigma_g)
                except ValueError as err:
                    print(f"Could not converge for density:\nrho_ice = {ice_density}\nrho_sil = {sil_density}")
                    print(err)

                time.sleep(1)
                print(f"{' radius (cm) ':=^20}{' mass (g) ':=^20}{' density (g/cm3) ':=^20}{' ice fraction (%) ':=^20}{' W(a)da (g/cm2) ':=^20}{' f(a)da (g/cm3) ':=^20}")
                for i in range(pebbles.n_pebbles):
                    print(f"{pebbles.radius[i]:^20.3e}{pebbles.mass[i]:^20.3e}{pebbles.density[i]:^20.5f}{pebbles.ice_fraction[i] * 100:^20.5f}{W_da[i]:^20.5f}{f_da[i]:^20.3e}")
                print("")

