import numpy as np

t = 'time'
p = 'proton'
a = 'alpha'

def scalar_velocity(
        data,
):

    data[p]['v_mag'] = []
    data[a]['v_mag'] = []

    L = len(data[p][t])

    for i in range(L):
        arg_p_ = (data[p]['vp1_x'][i]) ** 2 + (data[p]['vp1_y'][i]) ** 2 + (data[p]['vp1_z'][i]) ** 2
        arg_a_ = (data[a]['va_x'][i]) ** 2 + (data[a]['va_y'][i]) ** 2 + (data[a]['va_z'][i]) ** 2

        if arg_p_ < 0:
            arg_p_ = 0
        if arg_a_ < 0:
            arg_a_ = 0

        data[p]['v_mag'].append(np.sqrt(arg_p_))
        data[a]['v_mag'].append(np.sqrt(arg_a_))

    return data



def validate_theta(
        theta,
):
    for i in range(len(theta)):
        if theta[i] > 15:
            theta[i] = 15
        elif theta[i] < 0:
            theta[i] = 0
        else:
            pass

    return theta



def scalar_temps(
        solar_data,
        spc_data,
):
    key_names = {}
    file_names = []
    factor = 11604
    length = []
    for file in solar_data:
        file_names.append(file)
        key_names[file] = []
        for key in file:
            key_names[file] = solar_data[file].keys()

    L = len(solar_data[p][t])
    psp_result = {}
    psp_result_keys = ['proton_scalar_temp_1', 'proton_scalar_temp_2', 'alpha_scalar_temp', 'theta_ap', 'dens_ap',
                   'proton_perpar', 'alpha_perpar', 'proton_1_k', 'alpha_k','time','proton_R', 'alpha_R']

    for res_key in psp_result_keys:
        psp_result[res_key] = np.zeros(L)

    for i in range(L):
        psp_result['proton_scalar_temp_1'][i] = ((2 * solar_data[p]['Tperp1'][i] + solar_data[p]['Trat1'][i]) / 3)
        psp_result['proton_scalar_temp_2'][i] = ((2 * solar_data[p]['Tperp2'][i] + solar_data[p]['Trat2'][i]) / 3)
        psp_result['proton_perpar'][i] = solar_data[p]['Tperp1'][i] / solar_data[p]['Trat1'][i]
        psp_result['proton_1_k'][i] = psp_result['proton_scalar_temp_1'][i] * factor
        psp_result['proton_R'][i] = solar_data[p]['Tperp1'][i]/solar_data[p]['Trat1'][i]

        psp_result['alpha_scalar_temp'][i] = (2 * solar_data[a]['Ta_perp'][i] + solar_data[a]['Trat'][i]) / 3
        psp_result['alpha_perpar'][i] = solar_data[a]['Ta_perp'][i] / solar_data[a]['Trat'][i]
        psp_result['alpha_k'][i] = psp_result['alpha_scalar_temp'][i] * factor
        psp_result['alpha_R'][i] =solar_data[a]['Ta_perp'][i]/solar_data[a]['Trat'][i]

        if psp_result['proton_scalar_temp_1'][i] == 0:
            psp_result['theta_ap'][i] = 0
        else:
            psp_result['theta_ap'][i] = psp_result['alpha_scalar_temp'][i] / psp_result['proton_scalar_temp_1'][i]

        if solar_data[p]['np1'][i] == 0:
            psp_result['dens_ap'][i] = 0
        else:
            psp_result['dens_ap'][i] = solar_data[a]['na'][i] / solar_data[p]['np1'][i]

    psp_result['theta_ap'] = validate_theta(psp_result['theta_ap'])


    wind_result_keys = ['wind_alpha_scalar_temp', 'wind_proton_scalar_temp', 'wind_theta']
    wind_result = {}

    for i in range(len(wind_result_keys)):
        wind_result[wind_result_keys[i]] = np.zeros(len(spc_data['Wind_Temps.csv']['time']))

    for i in range(len(spc_data['Wind_Temps.csv']['time'])):
        arg_ = spc_data['Wind_Temps.csv']['TEMP_ALPHA_S/C_eV'][i]

        if arg_ == 0:
            wind_result[wind_result_keys[0]][i] = 10**6
        else:
            wind_result[wind_result_keys[0]][i] = arg_

        wind_result[wind_result_keys[1]][i] = spc_data['Wind_Temps.csv']['TEMP_PROTN_S/C_eV'][i]

        wind_result['wind_theta'][i] = spc_data['Wind_Temps.csv']['TEMP_ALPHA_S/C_eV'][i]/spc_data['Wind_Temps.csv']['TEMP_PROTN_S/C_eV'][i]

    wind_result['wind_theta'] = validate_theta(wind_result['wind_theta'])


    return psp_result, wind_result



