import numpy as np
import random
from modules.collisions.model import theta_ap as theta_ap_
from modules.core.constants import const as core_const

t = 'time'
p = 'proton'
a = 'alpha'


def scalar_velocity(
        data,
):
    for key in data.keys():
        encount = key
        file_val = []

        for value in data[key].keys():
            file_val.append(value)

        p = file_val[0]
        a = file_val[1]

        data[encount][p]['v_mag'] = []
        data[encount][a]['v_mag'] = []

        L = len(data[encount][p][t])
        l = len(data[encount][a][t])

        for i in range(L):
            arg_p_ = (data[encount][p]['vp1_x'][i]) ** 2 + (
            data[encount][p]['vp1_y'][i]) ** 2 + (
                         data[encount][p]['vp1_z'][i]) ** 2
            if arg_p_ < 0:
                arg_p_ = 0
            data[encount][p]['v_mag'].append(np.sqrt(arg_p_))
        for i in range(l):
            arg_a_ = (data[encount][a]['va_x'][i]) ** 2 + (
            data[encount][a]['va_y'][i]) ** 2 + (
                         data[encount][a]['va_z'][i]) ** 2
            if arg_a_ < 0:
                arg_a_ = 0
            data[encount][a]['v_mag'].append(np.sqrt(arg_a_))

    return data


def scalar_temps(
        solar_data,
        spc_data,
):

    psp_res = {}
    wind_res = {}

    for encounter in solar_data.keys():
        x, y = temp_generate(solar_data[encounter], spc_data[encounter])
        psp_res[encounter] = x
        wind_res[encounter] = y

    return psp_res, wind_res


def temp_generate(
        solar_encounter,
        spc_encounter,
):
    factor = 11604

    psp_result = {}
    psp_result_keys = ['proton_scalar_temp_1', 'proton_scalar_temp_2',
                       'alpha_scalar_temp', 'theta_ap', 'dens_ap',
                       'proton_perpar', 'alpha_perpar', 'proton_1_k', 'alpha_k', 'time',
                       'proton_R', 'alpha_R', ]
    wind_result = {}
    wind_result_keys = ['wind_alpha_scalar_temp', 'wind_proton_scalar_temp', 'wind_theta']

    file_val = []
    for file in solar_encounter.keys():
        file_val.append(file)
    p = file_val[0]
    a = file_val[1]

    L = len(solar_encounter[p][t])

    for res_key in psp_result_keys:
        psp_result[res_key] = np.zeros(L)

    for res_key in wind_result_keys:
        wind_result[res_key] = np.zeros(L)

    for i in range(L):

        # print(solar_data[encount][p]['Tperp1'][i], solar_data[encount][p]['Trat1'][i])
        psp_result['proton_scalar_temp_1'][i] = (
                (2 * solar_encounter[p]['Tperp1'][i] + solar_encounter[p]['Trat1'][
                    i]) / 3)
        psp_result['proton_scalar_temp_2'][i] = (
                (2 * solar_encounter[p]['Tperp2'][i] + solar_encounter[p]['Trat2'][
                    i]) / 3)
        psp_result['proton_perpar'][i] = (solar_encounter[p]['Tperp1'][i] / \
                                          solar_encounter[p]['Trat1'][i])
        psp_result['proton_1_k'][i] = (psp_result['proton_scalar_temp_1'][i] * factor)
        psp_result['proton_R'][i] = (
                    solar_encounter[p]['Tperp1'][i] / solar_encounter[p]['Trat1'][i])

        psp_result['alpha_scalar_temp'][i] = ((2 * solar_encounter[a]['Ta_perp'][i] +
                                               solar_encounter[a]['Trat'][i]) / 3)
        psp_result['alpha_perpar'][i] = (solar_encounter[a]['Ta_perp'][i] / \
                                         solar_encounter[a]['Trat'][i])
        psp_result['alpha_k'][i] = (psp_result['alpha_scalar_temp'][i] * factor)
        psp_result['alpha_R'][i] = (
                    solar_encounter[a]['Ta_perp'][i] / solar_encounter[a]['Trat'][i])

        if psp_result['proton_scalar_temp_1'][i] == 0:
            psp_result['theta_ap'][i] = 0
        else:
            alpha_noise = random.randint(0,15)/100
            proton_noise = random.randint(0,7)/100
            val = random.random()
            if val > 0.5:
                alpha_noise = alpha_noise + 1
                proton_noise = proton_noise + 1
            else:
                alpha_noise = 1 - alpha_noise
                proton_noise = 1 - proton_noise

            psp_result['theta_ap'][i] = (psp_result['alpha_scalar_temp'][i] / \
                                         psp_result['proton_scalar_temp_1'][i])

        if solar_encounter[p]['np1'][i] == 0:
            psp_result['dens_ap'][i] = 0
        else:
            psp_result['dens_ap'][i] = (
                        solar_encounter[a]['na'][i] / solar_encounter[p]['np1'][i])

        wind = spc_encounter['Wind_Temps.csv']

        arg_ = spc_encounter['Wind_Temps.csv']['TEMP_ALPHA_S/C_eV'][i]

        if arg_ == 0:
            wind_result['wind_alpha_scalar_temp'][i] = 10 ** 6
        else:
            wind_result['wind_alpha_scalar_temp'][i] = arg_

        wind_result['wind_proton_scalar_temp'][i] = (wind['TEMP_PROTN_S/C_eV'][i])

        if wind['TEMP_PROTN_S/C_eV'][i] == 0:
            wind['TEMP_PROTN_S/C_eV'][i] = 10 * 30

        wind_result['wind_theta'][i] = (
                    wind['TEMP_ALPHA_S/C_eV'][i] / wind['TEMP_PROTN_S/C_eV'][i])

        # psp_result['theta_ap'] = theta_ap_.validate_theta(psp_result['theta_ap'])
        # wind_result['wind_theta'] = theta_ap_.validate_theta(wind_result['wind_theta'])

        res_psp = psp_result
        res_wind = wind_result

    return res_psp, res_wind
