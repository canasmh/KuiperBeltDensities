from pebble_accretion import focus_accretion, bondi_accretion, hill_accretion


def rk4(dt, acc_type, chi=0.4, gamma=0.65, alpha=1e-4, **kwargs):
    if acc_type.upper() == "FOCUS":
        needed_keys = ["m", "rho", "rho_d", "Omega", "St", "deltav", "H_g"]
        for key in needed_keys:
            if key not in kwargs.keys():
                raise TypeError(f"missing required keyword argument '{key}'")
        
        m = kwargs[needed_keys[0]]
        rho = kwargs[needed_keys[1]]
        rho_d = kwargs[needed_keys[2]]
        Omega = kwargs[needed_keys[3]]
        St = kwargs[needed_keys[4]]
        deltav = kwargs[needed_keys[5]]
        H_g = kwargs[needed_keys[6]]

        m1 = focus_accretion(m, rho, rho_d, Omega, St, deltav, H_g, alpha)
        m2 = focus_accretion(m + 0.5 * dt * m1, rho, rho_d, Omega, St, deltav, H_g, alpha)
        m3 = focus_accretion(m + 0.5 * dt * m2, rho, rho_d, Omega, St, deltav, H_g, alpha)
        m4 = focus_accretion(m + dt * m3, rho, rho_d, Omega, St, deltav, H_g, alpha)
        dm = dt / 6 * (m1 + 2 * m2 + 2 * m3 + m4)
    
    elif acc_type.upper() == "BONDI":
        needed_keys = ["m", "r", "rho_d", "Omega", "St", "deltav", "H_g"]
        for key in needed_keys:
            if key not in kwargs.keys():
                raise TypeError(f"missing required keyword argument '{key}'")
        
        m = kwargs[needed_keys[0]]
        r = kwargs[needed_keys[1]]
        rho_d = kwargs[needed_keys[2]]
        Omega = kwargs[needed_keys[3]]
        St = kwargs[needed_keys[4]]
        deltav = kwargs[needed_keys[5]]
        H_g = kwargs[needed_keys[6]]

        m1 = bondi_accretion(m, r, rho_d, Omega, St, deltav, H_g, chi, gamma, alpha)
        m2 = bondi_accretion(m + 0.5 * dt * m1, r, rho_d, Omega, St, deltav, H_g, chi, gamma, alpha)
        m3 = bondi_accretion(m + 0.5 * dt * m2, r, rho_d, Omega, St, deltav, H_g, chi, gamma, alpha)
        m4 = bondi_accretion(m + dt * m3, r, rho_d, Omega, St, deltav, H_g, chi, gamma, alpha)
        dm = dt / 6 * (m1 + 2 * m2 + 2 * m3 + m4)

    elif acc_type.upper() == "HILL":
        needed_keys = ["m", "r", "rho_d", "Omega", "St", "deltav", "H_g"]
        for key in needed_keys:
            if key not in kwargs.keys():
                raise TypeError(f"missing required keyword argument '{key}'")
        
        m = kwargs[needed_keys[0]]
        r = kwargs[needed_keys[1]]
        rho_d = kwargs[needed_keys[2]]
        Omega = kwargs[needed_keys[3]]
        St = kwargs[needed_keys[4]]
        deltav = kwargs[needed_keys[5]]
        H_g = kwargs[needed_keys[6]]

        m1 = hill_accretion(m, r, rho_d, Omega, St, deltav, H_g, chi, gamma, alpha)
        m2 = hill_accretion(m + 0.5 * dt * m1, r, rho_d, Omega, St, deltav, H_g, chi, gamma, alpha)
        m3 = hill_accretion(m + 0.5 * dt * m2, r, rho_d, Omega, St, deltav, H_g, chi, gamma, alpha)
        m4 = hill_accretion(m + dt * m3, r, rho_d, Omega, St, deltav,  H_g, chi, gamma, alpha)
        dm = dt / 6 * (m1 + 2 * m2 + 2 * m3 + m4)

    else:
        raise TypeError("acc_type accepts 3 values: 'focus', 'bondi, and 'hill'")

    return dm