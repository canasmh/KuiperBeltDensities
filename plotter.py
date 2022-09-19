from cProfile import label
import matplotlib.pyplot as plt
import numpy as np

def plot_kbo_data(t):
    from reader import KuiperBeltData
    from constants import M_PLUTO

    SMALL = 10
    MEDIUM = 12
    LARGE = 14

    kbos = KuiperBeltData()

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)
    ax.set_title(f"T = {t :.4e} Years", fontsize=LARGE)

    ax.scatter(kbos.mass / M_PLUTO, kbos.density, marker="*", s=15 ** 2, c="r", zorder=2.5)
    ax.errorbar(x=kbos.mass / M_PLUTO, y=kbos.density, yerr=[kbos.min_density, kbos.max_density], ls='none', ecolor='k')

    ax.set_xscale('log')
    ax.set_xlabel(r'Mass (M$_{\rm{Pluto}}$)', fontsize=MEDIUM, labelpad=SMALL)
    ax.set_ylabel(r'Density (g cm$^{-3}$)', fontsize=MEDIUM, labelpad=SMALL)
    ax.set_xlim(6e-5, 2.5e0)
    ax.set_ylim(0, 3)
    plt.show()


if __name__ == "__main__":

    plot_kbo_data(t=0)

