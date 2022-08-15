import numpy as np
import matplotlib.pyplot as plt

from modules.core.loadsave import file_dir as fd
from modules.core.variables import string_man as sm, char_man as cm

slash = fd.slash()
path = fd.dir_path()
file_type = '.txt'
font = 'Courier New'

def fig_5(

):
    all_txt = fd.all_file_type(file_type, path)
    file_names = []

    for entry in all_txt:
        arg_ = sm.split(entry, slash)[-1]
        for i in range(len(file_type)):
            arg_ = cm.remove_end(arg_)
        entry = arg_
        file_names.append(entry)

    res = {'radius': {}, 'theta': {}}

    for name in file_names:
        value = name[0]
        if value == 't':
            number = name[5:len(name)]
        elif value == 'r':
            number = name[6:len(name)]
        else:
            raise ValueError(
                ''
            )

        res['radius'][number] = np.loadtxt(sm.slash_check(path) + 'radius' + str(number) + file_type)
        res['theta'][number] = np.loadtxt(sm.slash_check(path) + 'theta' + str(number) + file_type)

    L = len(res['theta'])
    ls = ['-', '-', '-.', '-', '-.', '--', '-.', ':', ':', '--', '-.', ':',]
    plt.figure(figsize=(10,10))
    j = 0
    for i in res['theta']:
        label = round(float(i), 4)
        plt.plot(res['radius'][i], res['theta'][i], linestyle=ls[j], linewidth=1, label=label)
        j=j+1

    plt.tight_layout(pad=10 * 0.8)
    plt.xticks(fontsize=24)
    plt.yticks(fontsize=24)
    plt.legend(loc='upper right')
    plt.grid()
    plt.ylabel(r'$\theta_{\alpha p} = \frac{T_{\alpha}}{T_{p}}$', fontsize=30, fontname=font)
    plt.xlabel(r'Radius ${\rm au}$', fontsize=30, fontname=font)
    plt.show()

    return

fig_5()

