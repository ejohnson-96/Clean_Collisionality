import scipy.stats

from modules.core.constants import const
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors, markers as mark
import random
import numpy as np
import warnings
from modules.core.variables import char_man as cm
from modules.core.features import smooth as smoothing
from modules.core.constants import const

const()


def colours_validate(
        colours,
):
    if colours is None:
        raise ValueError(
            "Error: No colours were provided."
        )
    else:

        if not isinstance(colours, list):
            raise TypeError(
                "Error: Colours must be a list, instead "
                f"got type {type(colours)}."
            )

        colors = dict(mcolors.BASE_COLORS, **mcolors.CSS4_COLORS)
        by_hsv = sorted((tuple(mcolors.rgb_to_hsv(mcolors.to_rgba(color)[:3])), name)
                        for name, color in colors.items())
        valid_colors = [name for hsv, name in by_hsv]
        res = []

        for color in colours:
            if not isinstance(color, str):
                raise ValueError(
                    f"Error: The entry {color} in colours, is not a "
                    f"string, instead got type {type(color)}."
                )
            elif color not in valid_colors:
                raise ValueError(
                    f"Error: The entry {color} in colours is not "
                    "supported, please try another colour."
                )
            else:
                res.append(color)

    return res


def validate_styles(
        line_style,
):
    valid_styles = []
    for style in mark.MarkerStyle.markers.keys():
        valid_styles.append(style)

    if isinstance(line_style, str):
        if line_style in valid_styles:
            return line_style
        else:
            raise ValueError(
                f"Error: Line style time {line_style} is not supported."
            )
    elif isinstance(line_style, list):
        for style in line_style:
            if not isinstance(style, str):
                raise TypeError(
                    f"Error: Style '{style}' in line styles is not a "
                    f"string, instead got type {type(style)}."

                )
            else:
                return line_style
    else:
        raise TypeError(
            f"Error: Line styles must be a string or a list of strings,"
            f" instead got type of {type(line_style)}."
        )


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
        colours=None,
        style_line='-',

):
    plt.figure(figsize=(const.x_dim, const.y_dim))

    if not isinstance(x_data, (list, np.ndarray)):
        raise ValueError(
            f"Error: x_data type not supported, {type(x_data)}"
        )
    x_0 = x_data[0]
    y_0 = 0

    style_line = validate_styles(style_line)
    if isinstance(y_data, (list, np.ndarray)):
        plt.plot(x_data, y_data, label=cm.capital_first_letter(label), color='black',
                 linewidth=const.line_width, linestyle=style_line)
    elif isinstance(y_data, dict):
        y_labels = []
        for key in y_data.keys():
            y_labels.append(key)
        for i in range(len(y_labels)):
            if colours is None:
                r = random.random()
                g = random.random()
                b = random.random()
                if isinstance(style_line, str):
                    plt.plot(x_data, y_data[y_labels[i]],
                             label=cm.capital_first_letter(y_labels[i]), color=(r, g, b),
                             linewidth=const.line_width, linestyle=style_line)
                else:
                    plt.plot(x_data, y_data[y_labels[i]],
                             label=cm.capital_first_letter(y_labels[i]), color=(r, g, b),
                             linewidth=const.line_width, linestyle=style_line[i])
            else:
                if not isinstance(colours, (list, np.ndarray)):
                    raise TypeError(
                        "Error: Colours needs to be list or an array,"
                        f"instead, got type {type(colours)}."
                    )
                elif len(colours) != len(y_data):
                    raise ValueError(
                        "Error: Colours must have the same length as the"
                        "number of y data entries, instead got lengths "
                        f"of {len(colours)} and {len(y_data)} respectively."
                    )
                else:
                    colours = colours_validate(colours)
                    colour = colours[i]
                    if isinstance(style_line, str):
                        plt.plot(x_data, y_data[y_labels[i]],
                                 label=cm.capital_first_letter(y_labels[i]), color=colour,
                                 linewidth=const.line_width, linestyle=style_line)

                    else:
                        plt.plot(x_data, y_data[y_labels[i]],
                                 label=cm.capital_first_letter(y_labels[i]), color=colour,
                                 linewidth=const.line_width, linestyle=style_line[i])

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
        colours=None,
        style='-',
):
    if isinstance(y_data, dict):
        y_arg_ = {}
        y_ = {}
        for key in y_data:
            y_arg_[key] = smoothing.smooth(y_data[key])
            bin_num = int(max(y_arg_[key]) - min(y_arg_[key]) / const.bin_width)

            hist = np.histogram(y_arg_[key], bins=bin_num)
            hist_dist = scipy.stats.rv_histogram(hist)
            y_[key] = smoothing.smooth(hist_dist.pdf(x_data), const.pdf_smooth)
    elif isinstance(y_data, (list, np.ndarray)):
        y_arg_ = smoothing.smooth(y_data)
        bin_num = int((max(y_arg_) - min(y_arg_)) / const.bin_width)
        hist = np.histogram(y_arg_, bins=bin_num)
        hist_dist = scipy.stats.rv_histogram(hist)
        y_ = smoothing.smooth(hist_dist.pdf(x_data), const.pdf_smooth)

    else:
        raise TypeError(
            "Error: Argument must be a dict, list or array instead "
            f"got type {type(y_data)}."
        )

    graph(x_data, y_, x_lim, y_lim, limits, degree, title, label, x_axis, y_axis, grid,
          y_log, x_log, colours, style)

    return
