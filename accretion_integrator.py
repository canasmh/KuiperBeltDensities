

from pebble_accretion import focus_accretion


def rk4(dt, acc_type, alpha=1e-4, **kwargs):
    if acc_type.upper() == "FOCUS":
        needed_keys = ["m", "rho", "rho_d", "Omega", "St", "vesc", "deltav", "H_g"]
        for key in needed_keys:
            if key not in kwargs.keys():
                raise TypeError(f"missing required keyword argument '{key}'")
        
        m = kwargs[needed_keys[0]]
        rho = kwargs[needed_keys[1]]
        rho_d = kwargs[needed_keys[2]]
        Omega = kwargs[needed_keys[3]]
        St = kwargs[needed_keys[4]]
        vesc = kwargs[needed_keys[5]]
        deltav = kwargs[needed_keys[6]]
        H_g = kwargs[needed_keys[7]]

        m1 = focus_accretion(m, rho, rho_d, Omega, St, vesc, deltav, H_g, alpha)
        m2 = focus_accretion(m + 0.5 * dt * m1, rho, rho_d, Omega, St, vesc, deltav, H_g, alpha)
        m3 = focus_accretion(m + 0.5 * dt * m2, rho, rho_d, Omega, St, vesc, deltav, H_g, alpha)
        m4 = focus_accretion(m + dt * m3, rho, rho_d, Omega, St, vesc, deltav, H_g, alpha)
        dm = dt / 6 * (m1 + 2 * m2 + 2 * m3 + m4)
    
    elif acc_type.upper() == "BONDI":
        pass

    elif acc_type.upper() == "HILL":
        pass

    else:
        raise ("Need to declare acc type")

    return dm