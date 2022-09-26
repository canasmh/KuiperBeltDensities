from cProfile import label
import matplotlib.pyplot as plt
from reader import StreamingInstabilityData
import numpy as np

def plot_kbo_data(t, si_data : StreamingInstabilityData, i=None, savefig=True):
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
    img = ax.scatter(x=si_data.mass / M_PLUTO, y=si_data.density, c=si_data.ice_fraction * 100, vmin=0, vmax=100)

    color_bars = np.linspace(0, 100, 11)
    cb_label = make_string(color_bars)

    cbar = fig.colorbar(img, ax=ax, extendrect=True, drawedges=False,
                        ticks=color_bars)
    cbar.set_label("Ice Fraction (%)", fontsize=12, labelpad=20, rotation=270)
    cbar.ax.set_yticklabels(cb_label, fontsize=12)

    ax.set_xscale('log')
    ax.set_xlabel(r'Mass (M$_{\rm{Pluto}}$)', fontsize=MEDIUM, labelpad=SMALL)
    ax.set_ylabel(r'Density (g cm$^{-3}$)', fontsize=MEDIUM, labelpad=SMALL)
    ax.set_xlim(6e-5, 2.5e0)
    ax.set_ylim(0, 3)

    if savefig:
        plt.savefig('./images/_tmp%04d.png' % i, dpi=350)
        plt.close()
    else:
        plt.show()


def make_string(some_list):

    return [f"{item:.1f}" for item in some_list]

if __name__ == "__main__":

    kbos = StreamingInstabilityData(rho_ice=1.0, rho_sil=3.0, unit_mass=2.823973078884959e+28)
    plot_kbo_data(t=0, si_data=kbos, savefig=False)

