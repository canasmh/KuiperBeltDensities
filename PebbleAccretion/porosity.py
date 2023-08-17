"""
    The following porosity modeled was found be applying a 6th order polynomial fit
    to the BN19 models:

    For f_m = 0.5:
        y = 5.71165285e-18 * x ** 6
          - 2.11541836e-14 * x ** 5
          + 2.71681541e-11 * x ** 4 
          - 1.31982250e-08 * x ** 3 
          + 1.03115299e-06 * x ** 2
          - 6.27393765e-05 * x ** 1
          + 6.00049944e-01
    
    For f_m = 0.7:
        y = 2.06933274e-17 * x ** 6
          - 7.21749477e-14 * x ** 5
          + 8.79474892e-11 * x ** 4 
          - 4.33002256e-08 * x ** 3 
          + 6.76651357e-06 * x ** 2
          - 5.49232826e-04 * x ** 1
          + 6.12671290e-01

    For f_m = 0.8:
        y = 2.09194677e-18 * x ** 6
          + 2.48659476e-15 * x ** 5
          - 2.45239147e-11 * x ** 4 
          + 3.57165654e-08 * x ** 3 
          - 1.88578106e-05 * x ** 2
          + 2.57076850e-03 * x ** 1
          + 5.05481730e-01

    But do the polynomial nature of the graph, need to apply some additional conditions.
    Apply a constant below some diameter:
    for f_m = 0.5:
        if d <= 100:
            porosity = 0.593
        else if d >=734:
            porosity = -1.52908068e-04 * x - 0.05063
    for f_m = 0.7:
        if d <= 100:
            porosity = 0.5898
        else if d >=709:
            porosity = -9.48782536e-05 * x - 0.04573
    for f_m = 0.8:
        if d <= 140:
            porosity = 0.586
        else if d >=690:
            porosity = -6.71755725e-05 * x - 0.04165
"""


def get_porosity(radius, ice_fraction):
    if ice_fraction <= 0.5:
        if radius <= 100:
            porosity = 0.593
        elif radius >= 734:
            porosity = -1.52908068e-04 * radius + 0.27503
        else:
            porosity = 5.71165285e-18 * radius ** 6 \
            - 2.11541836e-14 * radius ** 5 \
            + 2.71681541e-11 * radius ** 4 \
            - 1.31982250e-08 * radius ** 3 \
            + 1.03115299e-06 * radius ** 2 \
            - 6.27393765e-05 * radius ** 1 \
            + 6.00049944e-01 
    elif ice_fraction <= 0.7:
        if radius <= 100:
            porosity = 0.5898
        elif radius >= 709:
            porosity = -9.48782536e-05 * radius + 0.18027
        else:
            porosity = 2.06933274e-17 * radius ** 6 \
            - 7.21749477e-14 * radius ** 5 \
            + 8.79474892e-11 * radius ** 4 \
            - 4.33002256e-08 * radius ** 3 \
            + 6.76651357e-06 * radius ** 2 \
            - 5.49232826e-04 * radius ** 1 \
            + 6.12671290e-01
    else:
        if radius <= 140:
            porosity = 0.586
        elif radius >= 609:
            porosity = -6.71755725e-05 * radius + 0.13435
        else:
            porosity = 2.09194677e-18 * radius ** 6 \
            + 2.48659476e-15 * radius ** 5 \
            - 2.45239147e-11 * radius ** 4 \
            + 3.57165654e-08 * radius ** 3 \
            - 1.88578106e-05 * radius ** 2 \
            + 2.57076850e-03 * radius ** 1 \
            + 5.05481730e-01
    
    if porosity < 0:
        porosity = 0
            
    return porosity


if __name__ == "__main__":
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    # The diameter and densities pulled from BN19
    bn50 = pd.read_csv("./data/bn50f.csv")
    bn70 = pd.read_csv("./data/bn70f.csv")
    bn80 = pd.read_csv("./data/bn80f.csv")

    bn50diam = bn50['x']
    bn50dens = bn50[' y']
    bn70diam = bn70['x']
    bn70dens = bn70[' y']
    bn80diam = bn80['x']
    bn80dens = bn80[' y']

    # fully compact densities
    bn50max = 1456.3696783887217
    bn70max = 1898.3329598916055
    bn80max = 2242.108514418066
    # Calculate porosities
    bn50por = 1 - bn50dens / bn50max
    bn70por = 1 - bn70dens / bn70max
    bn80por = 1 - bn80dens / bn80max

    xnew = np.logspace(np.log10(20), np.log10(3000), 500)
    mc50por = []
    mc70por =[]
    mc80por = []

    for i in range(len(xnew)):
        mc50por.append(get_porosity(xnew[i] / 2, ice_fraction = 0.5))
        mc70por.append(get_porosity(xnew[i] / 2, ice_fraction = 0.7))
        mc80por.append(get_porosity(xnew[i] / 2, ice_fraction = 0.8))

    plt.figure(figsize=(10, 7))
    plt.scatter(bn50diam, bn50por, color='purple')
    plt.scatter(bn70diam, bn70por, color='green')
    plt.scatter(bn80diam, bn80por, color='orange')
    plt.plot(xnew, mc50por, c='purple', linestyle='--')
    plt.plot(xnew, mc70por, c='green', linestyle='--')
    plt.plot(xnew, mc80por, c='orange', linestyle='--')
    plt.xscale('log')
    plt.ylim(-0.05, 1.05)
    plt.plot()
    plt.show()
