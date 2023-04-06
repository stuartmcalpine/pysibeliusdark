import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

"""
Plot Figure 4 from Sibelius-DARK paper.

Plots mean density of dark matter in shells and spheres centred on the MW.

rho_mean used in computation = 127466270036.20662 * 0.307
"""

f,axarr = plt.subplots(2,1,figsize=(6,6))

for i in range(2):
    if i == 0:
        x,y,_ = np.loadtxt(f'./data/sibelius_dark_1.txt', unpack=True)
    else:
        x,y,_ = np.loadtxt(f'./data/sibelius_dark_1_cumulative.txt',unpack=True)
    axarr[i].plot(x, y, label='Sibelius-DARK',
            path_effects=[pe.Stroke(linewidth=2.5, foreground='k'), pe.Normal()],zorder=10)

# Some known cluster positions in Sibelius-DARK.
for name, dist, offset in zip(['Virgo', 'Perseus', 'Coma', 'Norma', 'Hercules'],
        [21.2, 77.2, 108.2,73.9,141.5],
        [0,0,0,-10,0]):
    axarr[0].arrow(dist, -1, 0., 0.15, head_width=2, head_length=0.075)
    axarr[0].text(dist+3+offset, -1 + 0.25, name, rotation=45, fontsize=6)

axarr[0].set_ylim(-1,1)
axarr[0].set_xlim(0,200)
axarr[1].set_xlim(0,500)
axarr[1].set_ylim(-0.5, 0.5)
axarr[0].set_ylabel(r'$(\rho(r) - \rho_{\mathrm{mean}}) / \rho_{\mathrm{mean}}$')
axarr[1].set_ylabel(r'$(\rho(\leq r) - \rho_{\mathrm{mean}}) / \rho_{\mathrm{mean}}$')
axarr[1].axvline(200, lw=0.5, ls='--', c='k')
axarr[1].legend(loc='upper right')

for i in range(2):
    axarr[i].axhline(0, lw=0.5, ls='--', c='k')

    axarr[i].set_xlabel(r'$d_{\mathrm{MW}}$ [Mpc]')
    axarr[i].yaxis.set_ticks_position('both')
    axarr[i].xaxis.set_ticks_position('both')
    axarr[i].minorticks_on()
plt.tight_layout(pad=0.1)
plt.savefig('density_vs_rad.pdf')
plt.close()

