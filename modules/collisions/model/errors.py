import numpy as np

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
    print(error_data.keys())
    for key in error_data.keys():
        for value in error_data[key].keys():
            print(value)
    proton_temp = scalar_temps['proton_scalar_temp_1']
    alpha_temp = scalar_temps['alpha_scalar_temp']
    p_s_e_x = error_data[p]['dV1x']
    p_s_e_y = error_data[p]['dV1y']
    p_s_e_z = error_data[p]['dV1z']
    a_s_e_x = error_data[a]['dVx']
    a_s_e_y = error_data[a]['dVy']
    a_s_e_z = error_data[a]['dVz']
    proton_speed = solar_data[p]['vmag']
    alpha_speed = solar_data[a]['vmag']
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
    proton_error_tpar = proton_error_isotropy * proton_error_tperp
    alpha_error_tpar = alpha_error_isotropy * alpha_error_tperp

    #par/pere

    proton_temp_error = np.sqrt((2*proton_error_tperp)**2+(proton_error_tpar**2))/3
    alpha_temp_error = np.sqrt((alpha_error_tperp)**2+(alpha_error_tpar)**2)/3

    proton_error_over_temp = proton_temp_error/proton_temp
    alpha_error_over_temp = alpha_temp_error/alpha_temp

    proton_error_over_dense = proton_error_density/proton_density
    alpha_error_over_dens = alpha_error_density/alpha_density

    proton_error_over_isotropy = proton_error_isotropy/proton_isotropy
    alpha_error_over_isotropy = alpha_error_isotropy/alpha_isotropy

    proton_error_speed = np.sqrt(p_s_e_x**2+p_s_e_y**2+p_s_e_z**2)
    alpha_error_speed = np.sqrt(a_s_e_z**2 + a_s_e_y**2 + a_s_e_x**2)

    proton_error_over_speed = proton_error_speed/proton_speed
    alpha_error_over_speed = alpha_error_speed/alpha_speed

    res = {}
    res[p] = {'Scalar Temp': proton_error_over_temp, 'Density': proton_error_over_dense, 'Speed': proton_error_speed, 'Isotropy': proton_error_over_isotropy}
    res[a] = {'Scalar Temp': alpha_error_over_temp, 'Density': alpha_error_over_dens, 'Speed': alpha_error_speed, 'Isotropy': alpha_error_over_isotropy}


    return res

#spacing for invtervals