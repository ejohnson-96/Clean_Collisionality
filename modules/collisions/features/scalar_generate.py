import numpy as np
from modules.collisions.model import theta_ap as theta_ap_
from modules.core.constants import const as core_const

t = 'time'
p = 'proton'
a = 'alpha'


def scalar_velocity(
        data,
):

    file_val = []
    for key in data.keys():
        encount = key
        for value in data[key].keys():
            file_val.append(value)
        p = file_val[0]
        a = file_val[1]

        data[encount][p]['v_mag'] = []
        data[encount][a]['v_mag'] = []

        L = len(data[encount][p][t])

        for i in range(L):
            arg_p_ = (data[encount][p]['vp1_x'][i]) ** 2 + (data[encount][p]['vp1_y'][i]) ** 2 + (
            data[encount][p]['vp1_z'][i]) ** 2
            arg_a_ = (data[encount][a]['va_x'][i]) ** 2 + (data[encount][a]['va_y'][i]) ** 2 + (
            data[encount][a]['va_z'][i]) ** 2

            if arg_p_ < 0:
                arg_p_ = 0
            if arg_a_ < 0:
                arg_a_ = 0

            data[encount][p]['v_mag'].append(np.sqrt(arg_p_))
            data[encount][a]['v_mag'].append(np.sqrt(arg_a_))

    return data


def scalar_temps(
        solar_data,
        spc_data,
):

    factor = 11604

    for key in solar_data.keys():
        for value in solar_data[key].keys():
            L = len(solar_data[key][value])

    psp_result = {}
    psp_result_keys = ['proton_scalar_temp_1', 'proton_scalar_temp_2',
                       'alpha_scalar_temp', 'theta_ap', 'dens_ap',
                       'proton_perpar', 'alpha_perpar', 'proton_1_k', 'alpha_k', 'time',
                       'proton_R', 'alpha_R',]



    wind_result_keys = ['wind_alpha_scalar_temp', 'wind_proton_scalar_temp', 'wind_theta']
    wind_result = {}

    res_psp = {}
    res_wind = {}

    for res_key in psp_result_keys:
        psp_result[res_key] = np.zeros(L)

    for i in range(len(wind_result_keys)):
        wind_result[wind_result_keys[i]] = np.zeros(L)

    file_val = []
    for key in solar_data.keys():
        encount = key
        wind = spc_data[encount]['Wind_Temps.csv']
        for value in solar_data[key].keys():
            file_val.append(value)
        p = file_val[0]
        a = file_val[1]

        for i in range(L):
            psp_result['proton_scalar_temp_1'][i] = (
                        (2 * solar_data[encount][p]['Tperp1'][i] + solar_data[encount][p]['Trat1'][i]) / 3)
            psp_result['proton_scalar_temp_2'][i] = (
                        (2 * solar_data[encount][p]['Tperp2'][i] + solar_data[encount][p]['Trat2'][i]) / 3)
            psp_result['proton_perpar'][i] = solar_data[encount][p]['Tperp1'][i] / \
                                             solar_data[encount][p]['Trat1'][i]
            psp_result['proton_1_k'][i] = psp_result['proton_scalar_temp_1'][i] * factor
            psp_result['proton_R'][i] = solar_data[encount][p]['Tperp1'][i] / solar_data[encount][p]['Trat1'][i]

            psp_result['alpha_scalar_temp'][i] = (2 * solar_data[encount][a]['Ta_perp'][i] +
                                                  solar_data[encount][a]['Trat'][i]) / 3
            psp_result['alpha_perpar'][i] = solar_data[encount][a]['Ta_perp'][i] / \
                                            solar_data[encount][a]['Trat'][i]
            psp_result['alpha_k'][i] = psp_result['alpha_scalar_temp'][i] * factor
            psp_result['alpha_R'][i] = solar_data[encount][a]['Ta_perp'][i] / solar_data[encount][a]['Trat'][i]

            if psp_result['proton_scalar_temp_1'][i] == 0:
                psp_result['theta_ap'][i] = 0
            else:
                psp_result['theta_ap'][i] = psp_result['alpha_scalar_temp'][i] / \
                                            psp_result['proton_scalar_temp_1'][i]

            if solar_data[encount][p]['np1'][i] == 0:
                psp_result['dens_ap'][i] = 0
            else:
                psp_result['dens_ap'][i] = solar_data[encount][a]['na'][i] / solar_data[encount][p]['np1'][i]
            psp_result['theta_ap'] = theta_ap_.validate_theta(psp_result['theta_ap'])
            res_psp[encount] = psp_result


            arg_ = spc_data[encount]['Wind_Temps.csv']['TEMP_ALPHA_S/C_eV'][i]

            if arg_ == 0:
                wind_result['wind_alpha_scalar_temp'][i] = 10 ** 6
            else:
                wind_result['wind_alpha_scalar_temp'][i] = arg_

            wind_result['wind_proton_scalar_temp'][i] = wind['TEMP_PROTN_S/C_eV'][i]

            if wind['TEMP_PROTN_S/C_eV'][i] == 0:
                wind['TEMP_PROTN_S/C_eV'][i] = 10*30

            wind_result['wind_theta'][i] = wind['TEMP_ALPHA_S/C_eV'][i] /wind['TEMP_PROTN_S/C_eV'][i]
            wind_result['wind_theta'] = theta_ap_.validate_theta(wind_result['wind_theta'])
            res_wind[encount] = wind_result


    return res_psp, res_wind
