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
        error_files,
):
    l = -1
    k = -1
    for particle in solar_data.keys():
        l = l + 1
        for key in solar_data[particle].keys():
            k = k + 1
            val = const.data_units[l][k]
            solar_data[particle][key] = clean(solar_data[particle][key], const.var_min[val],
                                         const.var_max[val])
        print(f"{(l / len(solar_data)) * 100:.2f} %", end="\r")
        k = -1


    if error_files:
        l = -1
        k = -1
        for particle in error_data.keys():
            l = l + 1
            for key in error_data[particle].keys():
                k = k + 1
                val = const.error_units[l][k]
                error_data[particle][key] = clean(error_data[particle][key], const.var_min[val],
                                             const.var_max[val])
            print(f"{(l / len(error_data)) * 100:.2f} %", end="\r")
            k = -1

    m = - 1
    n = - 1
    for orbit in spc_data.keys():
        m = m + 1
        for key in spc_data[orbit].keys():
            n = n + 1
            val = const.sc_units[m][n]
            spc_data[orbit][key] = clean(spc_data[orbit][key], const.var_min[val], const.var_max[val])
        print(f"{(m / len(spc_data)) * 100:.2f} %", end="\r")
        n = -1


    print('Data Scrub Complete', '\n')

    return solar_data, error_data, spc_data



def error_clean(
        solar_temps,
        sigma_temps,
):
    for i in range(len(solar_temps)):
        arg_ = sigma_temps[i]/solar_temps[i]


    return