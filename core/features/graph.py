from core.constants import const
import matplotlib.pyplot as plt
import random
import numpy as np
import warnings
from core.variables import char_man as cm


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

):
    plt.figure(figsize=(const.x_dim, const.y_dim))

    if not isinstance(x_data, (list, np.ndarray)):
        raise ValueError(
            f"Error: x_data type not supported, {type(x_data)}"
        )
    x_0 = x_data[0]
    y_0 = 0
    if isinstance(y_data, (list, np.ndarray)):
        plt.plot(x_data, y_data, label=cm.capital_first_letter(label), color='black', linewidth=const.line_width)
    elif isinstance(y_data, dict):
        y_labels = []
        for key in y_data.keys():
            y_labels.append(key)
        for i in range(len(y_labels)):
            r = random.random()
            g = random.random()
            b = random.random()
            plt.plot(x_data, y_data[y_labels[i]], label=cm.capital_first_letter(y_labels[i]), color=(r, g, b),
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

    arg_ = int(((const.x_dim + const.y_dim)/2))
    plt.tight_layout(pad=arg_*0.5)

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
    plt.xticks(rotation=degree)
    plt.show()

    return
