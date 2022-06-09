import scipy.stats

from modules.core.constants import const
import matplotlib.pyplot as plt
import random
import numpy as np
import warnings
from modules.core.variables import char_man as cm
from modules.core.features import smooth as smooth
from modules.core.constants import const

const()


def graph(
        x_data,
        y_data,
        x_lim=0,
        y_lim=0,
        limits=False,
        degree=0,
        title='Graph',
        label='Label',
        x_axis='x Data',
        y_axis='y Data',
        grid=True,
        y_log=False,
        x_log=False,

):
    plt.figure(figsize=(const.x_dim, const.y_dim))

    if not isinstance(x_data, (list, np.ndarray)):
        raise ValueError(
            f"Error: x_data type not supported, {type(x_data)}"
        )
    x_0 = x_data[0]
    y_0 = 0
    if isinstance(y_data, (list, np.ndarray)):
        plt.plot(x_data, y_data, label=cm.capital_first_letter(label), color='black',
                 linewidth=const.line_width)
    elif isinstance(y_data, dict):
        y_labels = []
        for key in y_data.keys():
            y_labels.append(key)
        for i in range(len(y_labels)):
            r = random.random()
            g = random.random()
            b = random.random()
            plt.plot(x_data, y_data[y_labels[i]],
                     label=cm.capital_first_letter(y_labels[i]), color=(r, g, b),
                     linewidth=const.line_width)
    else:
        raise ValueError(
            f"Error: y_data type not supported, {type(y_data)}"
        )

    plt.legend(loc='upper right', prop={'size': const.legend_size})
    plt.title(title,
              fontsize=const.title_size,
              fontname=const.font_family)
    plt.ylabel(y_axis, fontsize=const.label_size,
               fontname=const.font_family)
    plt.xlabel(x_axis, fontsize=const.label_size,
               fontname=const.font_family)
    plt.xticks(fontsize=const.tick_size)
    plt.yticks(fontsize=const.tick_size)

    arg_ = int(((const.x_dim + const.y_dim) / 2))
    plt.tight_layout(pad=arg_ * 0.5)

    if not limits:
        pass
    else:
        if x_0 > x_lim:
            warnings.warn(
                'x_0 is larger than x_lim, limits reverted to default values.'
                f"Values given; x_0: {x_0}, and x_lim: {x_lim}."
            )
        else:
            plt.xlim([x_0, x_lim])
        if y_0 > y_lim:
            warnings.warn(
                'y_0 is larger than y_lim, limits reverted to default values.'
                f"Values given; y_0: {y_0}, and y_lim: {y_lim}."
            )
        else:
            plt.ylim([y_0, y_lim])
    if not grid:
        pass
    else:
        plt.grid()
    if not x_log:
        pass
    else:
        plt.set_xscale('log')
    if not y_log:
        pass
    else:
        plt.set_yscale('log')

    plt.xticks(rotation=degree)
    plt.show()

    return


def histogram(
        x_data,
        y_data,
        x_lim=0,
        y_lim=0,
        limits=False,
        degree=0,
        title='Graph',
        label='Label',
        x_axis='x Data',
        y_axis='y Data',
        grid=True,
        y_log=False,
        x_log=False,

):
    y_arg_ = {}
    y_ = {}
    bins = []
    for key in y_data.key():
        y_arg_[key] = smooth(y_data[key])
        arg_ = int(max(y_arg_[key]) - min(y_arg_[key]) / const.bin_width)
        bins.append(arg_)

        hist = np.histogram(y_arg_[key])
        hist_dist = scipy.stats.rv_histogram(hist)

        y_[key] = smooth(hist_dist, const.pdf_smooth)

    graph(x_data, y_, x_lim, y_lim, limits, degree, title, label, x_axis, y_axis, grid,
          y_log, x_log, )

    return
