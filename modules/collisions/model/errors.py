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
