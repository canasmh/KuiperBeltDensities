from constants import AU_TO_CM, M_PLUTO, YRS_TO_SEC, G
import numpy as np
import matplotlib.pyplot as plt


from reader import KuiperBeltData
from constants import M_PLUTO

SMALL = 10
MEDIUM = 12
LARGE = 14

kbos = KuiperBeltData()

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
#ax.set_title(f"T = {t :.4e} Years", fontsize=LARGE)

ax.scatter(kbos.mass / M_PLUTO, kbos.density, marker="*", s=15 ** 2, c="r", zorder=2.5)
ax.errorbar(x=kbos.mass / M_PLUTO, y=kbos.density, yerr=[kbos.min_density, kbos.max_density], ls='none', ecolor='k')

xoffset=np.repeat( 0.   ,len(kbos.mass))
yoffset=np.repeat(-0.175,len(kbos.mass))

# Quaoar
xoffset[12]-=0.0375
yoffset[12]+=0.225

# Charon
xoffset[13] += 0.05
#yoffset[13]=

for i in range(len(kbos.mass)):
    name=str(kbos.name[i])
    if (name!='nan' and name!='Teharonhiawako'):
        ax.annotate(str(kbos.name[i]),((kbos.mass[i] / M_PLUTO)+xoffset[i], kbos.density[i]-kbos.min_density[i]+yoffset[i]),xycoords='data',ha='center',va='bottom')
        
#text, xy, xytext=None, xycoords='data',

ax.set_xscale('log')
ax.set_xlabel(r'Mass (M$_{\rm{Pluto}}$)', fontsize=MEDIUM, labelpad=SMALL)
ax.set_ylabel(r'Density (g cm$^{-3}$)', fontsize=MEDIUM, labelpad=SMALL)
ax.set_xlim(6e-5, 2.5e0)
ax.set_ylim(0, 3)

#plt.show()
plt.savefig('kbos_final.pdf')
