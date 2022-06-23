import numpy as np
from modules.collisions.constants import enc as encounter

t = 'time'
p = 'proton'
a = 'alpha'


def loop_uncer(
        solar_data,
        psp_scalar_temps,
):
    intervals_per_day = [i+1 for i in range(191)]
    interval_length = [24. / i for i in intervals_per_day]
    vals_ = {}
    particle_list = []
    for key in solar_data.keys():
        particle_list.append(key)
    particle_list.remove(t)

    for particle in particle_list:
        vals_[particle] = []

    for i in range(len(intervals_per_day)):
        arg_p_, arg_a_ = gen_uncer(solar_data, psp_scalar_temps, intervals_per_day[i])
        vals_[p].append(arg_p_)
        vals_[a].append(arg_a_)

    return interval_length, vals_


def gen_uncer(
        solar_data,
        scalar_temps,
        n_value,
):
    print(scalar_temps.keys(), solar_data.keys())
    L = len(solar_data[t])
    days = np.zeros(L)
    for i in range(L):
        days[i] = solar_data[p]['time'][i] / 86400  # 60

    temp_p = scalar_temps['proton_scalar_temp_1']  # solar_data[p]['v_mag']  #
    temp_a = scalar_temps['alpha_scalar_temp']  # solar_data[a]['v_mag'] #

    t_min = int(np.floor(min(days)))
    t_max = int(np.ceil(max(days)))
    n_days = t_max - t_min + 1
    n_interval = n_value

    out_mean_p = np.tile(0., [n_days, n_interval])
    out_mean_a = np.tile(0., [n_days, n_interval])

    for d in range(n_days):
        tk_d = np.where((days >= (d + t_min)) & (days < (d + t_min + 1)))[0]

        if len(tk_d) < 4:
            continue
        print(tk_d)
        tp_s = temp_p[tk_d]
        ta_s = temp_a[tk_d]
        dy_d = days[tk_d] - (d + t_min)

        for i in range(n_interval):
            tk_i = np.where((dy_d >= ((i + 0.) / n_interval)) &
                            (dy_d < ((i + 1.) / n_interval)))[0]
            # print(tk_i)
            if len(tk_i) < 2:
                continue

            out_mean_p[d, i] = np.std(tp_s[tk_i]) / np.mean(tp_s[tk_i])
            out_mean_a[d, i] = np.std(ta_s[tk_i]) / np.mean(ta_s[tk_i])
            # time[i] = d

    tk_p = np.where(out_mean_p > 0)
    tk_a = np.where(out_mean_a > 0)

    p_out = np.nanmedian(out_mean_p[tk_p]) * 100
    a_out = np.nanmedian(out_mean_a[tk_a]) * 100

    return p_out, a_out

def sigma_value(
        solar_data,
        error_data,
        scalar_temps,
):

    proton_temp = scalar_temps['proton_scalar_temp_1']
    alpha_temp = scalar_temps['alpha_scalar_temp']
    proton_temp_perp = solar_data[p]['Tperp1']
    alpha_temp_perp = solar_data[a]['Ta_perp']
    p_s_e_x = error_data[p]['dV1x']
    p_s_e_y = error_data[p]['dV1y']
    p_s_e_z = error_data[p]['dV1z']
    a_s_e_x = error_data[a]['dVx']
    a_s_e_y = error_data[a]['dVy']
    a_s_e_z = error_data[a]['dVz']
    proton_speed = solar_data[p]['v_mag']
    alpha_speed = solar_data[a]['v_mag']
    proton_density = solar_data[p]['np1']
    alpha_density = solar_data[a]['na']
    proton_isotropy = scalar_temps['proton_R']
    alpha_isotropy = scalar_temps['alpha_R']
    proton_error_tperp = error_data[p]['dT1_perp']
    alpha_error_tperp = error_data[a]['dT_perp']
    proton_error_density = error_data[p]['dn1']
    alpha_error_density = error_data[a]['dn']
    proton_error_isotropy = error_data[p]['dR1']
    alpha_error_isotropy = error_data[a]['dR']


    #par/pere
    L = len(solar_data['time'])
    proton_temp_error = np.zeros(L)
    alpha_temp_error = np.zeros(L)
    proton_error_over_temp = np.zeros(L)
    alpha_error_over_temp = np.zeros(L)
    proton_error_over_dense = np.zeros(L)
    alpha_error_over_dens = np.zeros(L)
    proton_error_over_isotropy = np.zeros(L)
    alpha_error_over_isotropy = np.zeros(L)
    proton_error_speed = np.zeros(L)
    alpha_error_speed = np.zeros(L)
    proton_error_over_speed = np.zeros(L)
    alpha_error_over_speed = np.zeros(L)
    proton_error_over_temp_perp = np.zeros(L)
    alpha_error_over_temp_perp = np.zeros(L)

    proton_error_tpar = np.zeros(L)
    alpha_error_tpar = np.zeros(L)

    fill_value = 10**31

    for i in range(L):
        proton_error_tpar[i] = proton_error_isotropy[i] * proton_error_tperp[i]
        alpha_error_tpar[i] = alpha_error_isotropy[i] * alpha_error_tperp[i]

        proton_temp_error[i] = np.sqrt((2*proton_error_tperp[i])**2+(proton_error_tpar[i]**2))/3
        alpha_temp_error[i] = np.sqrt((alpha_error_tperp[i])**2+(alpha_error_tpar[i])**2)/3

        if proton_temp[i] == 0:
            proton_temp[i] = fill_value
        proton_error_over_temp[i] = proton_temp_error[i]/proton_temp[i]
        if alpha_temp[i] == 0:
            alpha_temp[i] = fill_value
        alpha_error_over_temp[i] = alpha_temp_error[i]/alpha_temp[i]

        if proton_temp_perp[i] == 0:
            proton_temp_perp[i] = fill_value
        proton_error_over_temp_perp[i] = proton_error_tperp[i]/proton_temp_perp[i]
        if alpha_temp_perp[i] == 0:
            alpha_temp_perp[i] = fill_value
        alpha_error_over_temp_perp[i] = alpha_error_tperp[i]/alpha_temp_perp[i]

        if proton_density[i] == 0:
            proton_density[i] = fill_value
        proton_error_over_dense[i] = proton_error_density[i]/proton_density[i]
        if alpha_density[i] == 0:
            alpha_density[i] = fill_value
        alpha_error_over_dens[i] = alpha_error_density[i]/alpha_density[i]

        if proton_isotropy[i] == 0:
            proton_isotropy = fill_value
        proton_error_over_isotropy[i] = proton_error_isotropy[i]/proton_isotropy[i]
        if alpha_isotropy[i] == 0:
            alpha_isotropy[i] = fill_value
        alpha_error_over_isotropy[i] = alpha_error_isotropy[i]/alpha_isotropy[i]

        proton_error_speed[i] = np.sqrt(p_s_e_x[i]**2+p_s_e_y[i]**2+p_s_e_z[i]**2)
        alpha_error_speed[i] = np.sqrt(a_s_e_z[i]**2 + a_s_e_y[i]**2 + a_s_e_x[i]**2)

        proton_error_over_speed[i] = proton_error_speed[i]/proton_speed[i]
        alpha_error_over_speed[i] = alpha_error_speed[i]/alpha_speed[i]

    res = {}
    res[p] = {'Scalar Temp': proton_error_over_temp, 'Density': proton_error_over_dense, 'Speed': proton_error_over_speed, 'Isotropy': proton_error_over_isotropy, 'TPerp':proton_error_over_temp_perp}
    res[a] = {'Scalar Temp': alpha_error_over_temp, 'Density': alpha_error_over_dens, 'Speed': alpha_error_over_speed, 'Isotropy': alpha_error_over_isotropy, 'TPerp':alpha_error_over_temp_perp}


    return res

#spacing for invtervals