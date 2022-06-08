from modules.collisions.constants import *
from modules.collisions.features import lat_lon

const()


def clean(
        array,
        min_val=0,
        max_val=10 ** 30,
):
    for arg_name_ in (min_val, max_val):
        if not isinstance(arg_name_, (int, float)):
            raise TypeError(
                f"Argument '{arg_name_}' must be a number"
                f"instead get type '{type(arg_name_)}'"
            )

    for i in range(len(array)):
        array[i] = array[i]
        if abs(array[i]) > max_val:
            if i == 0:
                array[i] = 0
            else:
                array[i] = array[i - 1]
        elif abs(array[i]) < min_val:
            if i == 0:
                array[i] = 0
            else:
                array[i] = array[i - 1]
        else:
            array[i] = array[i]

    return array


def scrub_data(
        solar_data,
        error_data,
        spc_data,
):
    l = -1
    k = -1
    for x in solar_data.keys():
        l = l + 1
        for y in solar_data[x].keys():
            k = k + 1
            val = const.data_units[l][k]
            solar_data[x][y] = clean(solar_data[x][y], const.var_min[val],
                                     const.var_max[val])
        print(f"{(l / len(solar_data)) * 100:.2f} %", end="\r")
        k = -1
    m = -1
    n = -1

    l = -1
    k = -1
    for x in error_data.keys():
        l = l + 1
        for y in error_data[x].keys():
            k = k + 1
            val = const.error_units[l][k]
            error_data[x][y] = clean(error_data[x][y], const.var_min[val],
                                     const.var_max[val])
        print(f"{(l / len(error_data)) * 100:.2f} %", end="\r")
        k = -1
    m = -1
    n = -1

    for x in spc_data.keys():
        m = m + 1
        for y in spc_data[x].keys():
            n = n + 1
            val = const.sc_units[l][k]
            spc_data[x][y] = clean(spc_data[x][y], const.var_min[val], const.var_max[val])
        print(f"{(m / len(spc_data)) * 100:.2f} %", end="\r")
        n = -1
        m = -1

    print('Data Scrub Complete', '\n')

    return solar_data, error_data, spc_data