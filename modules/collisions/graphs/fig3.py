import numpy as np
import matplotlib.pyplot as plt

from modules.collisions.graphs.fig1 import gen_load_paths, make_bars
from modules.core.loadsave import file_dir as fd
from modules.core.variables import string_man as sm

paths = gen_load_paths()
slash = fd.slash()
font = 'Courier New'

def fig3(
        user_enc='EA',
        bin_number=100,
        value=0.07,
):

    path = paths[user_enc][1.0]

    theta10 = np.loadtxt(sm.jwos(path, 'theta_f.txt'))
    theta10 = theta10[np.logical_not(np.isnan(theta10))]
    thetaw = np.loadtxt(sm.jwos(path, 'wind_theta.txt'))
    thetaw = thetaw[np.logical_not(np.isnan(thetaw))]

    bn_ = bin_number
    val = value

    weg1 = np.ones_like(theta10) / float(len(theta10))
    weg2 = np.ones_like(thetaw) / float(len(thetaw))

    plt.figure(figsize=(10, 10))

    results, edges = np.histogram(theta10, range=(0, 14), weights=weg1, bins=bn_, density=1, )
    binWidth = edges[1] - edges[0]
    xs, ys = make_bars(edges[:-1], results * binWidth, val)
    plt.plot(xs, ys, color='black', linestyle='--', linewidth=1, label=r'1.0 ${\rm au}$', )

    results, edges = np.histogram(thetaw, range=(0, 14), weights=weg2, bins=bn_, density=1, )
    binWidth = edges[1] - edges[0]
    xs, ys = make_bars(edges[:-1], results * binWidth, val)
    plt.plot(xs, ys, color='blue', linestyle='-', linewidth=1, label=r'WIND @ 1.0 ${\rm au}$', )

    plt.tight_layout(pad=10 * 0.8)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.legend(loc='upper right', prop={'size': 26})
    plt.xlim(0, 14)
    plt.ylim(0, 0.05)
    plt.grid()
    plt.ylabel('Probability Density', fontsize=26, fontname=font)
    plt.xlabel(r'$\alpha$-Proton Relative Temperature, $\theta_{\alpha p} = \frac{T_{\alpha}}{T_{p}}$', fontsize=26,
               fontname=font)

    plt.show()


    return
