from cProfile import label
import matplotlib.pyplot as plt
from reader import StreamingInstabilityData
import numpy as np
from matplotlib import colors

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

    cmap = plt.cm.seismic_r
    norm = colors.BoundaryNorm(np.arange(0, 110, 10), cmap.N)

    ax.scatter(kbos.mass / M_PLUTO, kbos.density, marker="*", s=15 ** 2, c="darkgreen", zorder=2.5)
    ax.errorbar(x=kbos.mass / M_PLUTO, y=kbos.density, yerr=[kbos.min_density, kbos.max_density], ls='none', ecolor='darkgreen')
    img = ax.scatter(x=si_data.mass / M_PLUTO, y=si_data.density, c=si_data.ice_fraction * 100, cmap=cmap, norm=norm, zorder=3)

    color_bars = np.linspace(0, 100, 11)
    cb_label = make_string(color_bars)

    cbar = fig.colorbar(img, ax=ax, extendrect=True, drawedges=False,
                        ticks=color_bars)
    cbar.set_label("Ice Fraction (%)", fontsize=12, labelpad=20, rotation=270)
    cbar.ax.set_yticklabels(cb_label, fontsize=12)

    ax.set_xscale('log')
    ax.set_xlabel(r'Mass (M$_{\rm{Pluto}}$)', fontsize=MEDIUM, labelpad=SMALL)
    ax.set_ylabel(r'Density (g cm$^{-3}$)', fontsize=MEDIUM, labelpad=SMALL)
    ax.set_xlim(6e-5, 5)
    ax.set_ylim(0, 3)

    if savefig:
        plt.savefig('./images/_tmp%04d.png' % i, dpi=350, bbox_inches="tight")
        plt.close()
    else:
        plt.show()


def make_string(some_list):

    return [f"{item:.1f}" for item in some_list]


def animate_figures(filename, move_files=True):

    import os

    if move_files:
        dir_exists = os.system(f"mkdir images/{filename}")
        if dir_exists == 0:
            # This means directory did not exist already
            os.system(f"ffmpeg -i ./images/_tmp%04d.png -r 24 -b 1000k -vcodec mpeg4 -y ./images/{filename}/{filename}.mov")
            os.system(f"mv ./images/*.png ./images/{filename}")

        else:
            os.system(f"rm ./images/{filename}/*.png ./images/{filename}/{filename}.mov")
            os.system(f"ffmpeg -i ./images/_tmp%04d.png -r 24 -b 1000k -vcodec mpeg4 -y ./images/{filename}/{filename}.mov")
            os.system(f"mv ./images/*.png ./images/{filename}")
    

if __name__ == "__main__":

    kbos = StreamingInstabilityData(rho_ice=1.0, rho_sil=3.0, unit_mass=2.823973078884959e+28)
    plot_kbo_data(t=0, si_data=kbos, savefig=False)

