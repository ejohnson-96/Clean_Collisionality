import numpy as np
import matplotlib.pyplot as plt

from modules.core.loadsave import file_dir as fd
from modules.core.variables import string_man as sm, char_man as cm

slash = fd.slash()


def gen_load_paths(

):
    parent_name = fd.dir_name(fd.dir_parent())
    path = fd.dir_parent()
    dir_name = fd.dir_name(path)

    for i in range(len(parent_name) + len(dir_name) - 6):
        path = cm.remove_end(path)

    path = sm.jwos(sm.slash_check(path), 'data', slash, 'save', slash)

    encounters = ['EA', 'E4', 'E6', 'E7']
    vals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0]
    paths = {}
    for encount in encounters:
        paths[encount] = {}
        for val in vals:
            paths[encount][val] = sm.jwos(path, encount, slash, str(val), slash)

    return paths


paths = gen_load_paths()


def make_bars(
        x,
        y,
        width,
):
    xs = [x[0] - width]
    ys = [0]
    for i in range(len(x)):
        xs.append(x[i] - width)
        xs.append(x[i] + width)
        ys.append(y[i])
        ys.append(y[i])
    xs.append(x[-1] + width)
    ys.append(y[-1])

    return xs, ys


def fig1(
        user_enc='EA',
        bin_number=100,
        value=0.07,
):
    path = paths[user_enc][1.0]
    font = 'Courier New'

    theta = np.loadtxt(sm.jwos(path, 'theta_i.txt'))
    theta = theta[np.logical_not(np.isnan(theta))]

    bn_ = bin_number
    val = value
    weg2 = np.ones_like(theta) / float(len(theta))

    results, edges = np.histogram(theta, range=(0, 14), weights=weg2, bins=bn_, density=1, )
    binWidth = edges[1] - edges[0]
    xs, ys = make_bars(edges[:-1], results * binWidth, val)

    plt.figure(figsize=(10, 10))
    # plt.bar(edges[:-1], results*binWidth, binWidth)
    plt.plot(xs, ys, color='black', linestyle='--', linewidth=1, label=r'0.1 - 0.27 ${\rm au}$', )
    plt.tight_layout(pad=10 * 0.8)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.legend(loc='upper right', prop={'size': 32})
    plt.xlim(0, 14)
    plt.ylim(0, 0.05)
    plt.grid()
    plt.ylabel('Probability Density', fontsize=26, fontname=font)
    plt.xlabel(r'$\alpha$-Proton Relative Temperature, $\theta_{\alpha p} = \frac{T_{\alpha}}{T_{p}}$', fontsize=26, fontname=font)
    plt.show()
    plt.show()



    return

