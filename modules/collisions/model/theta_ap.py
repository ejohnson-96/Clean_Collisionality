import math
import numpy as np
from modules.collisions.constants import const as coll_const
from modules.core.constants import const as core_const
from modules.core.features import smooth as smoothing

coll_const()
core_const()

t = 'time'
p = 'proton'
a = 'alpha'
psp = 'PSP.csv'


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



def remove_theta(
        theta,
        value,
        rel_tol,
        smooth_=core_const.smooth,
):

    L = len(theta)
    if value is None:
        return theta
    else:
        if not isinstance(value, (float, int)):
            raise TypeError(
                "Error: Initial value must be a float or integer,"
                f"instead got type {type(value)}"
            )
        else:
            if not isinstance(rel_tol, (float, int)):
                raise TypeError(
                    "Error: Argument 'rel_tol' needs to be an integer"
                    f" or float, instead got type {type(rel_tol)}."
                )
            else:
                arg_min_ = value * (1 - rel_tol / 100)
                arg_max_ = value * (1 + rel_tol / 100)

                arg_ = theta[np.where((theta > arg_min_) & (theta < arg_max_))]
                try:
                    arg_v2 = [x for x in theta if x not in arg_]
                    arg_v3 = np.array(smoothing.smooth(arg_v2, smooth_))
                    res = np.resize(arg_v3, L)
                except:
                    res = np.array(smoothing.smooth(theta, smooth_))

                return res



def theta_ap_0(r_0, r_1, n_p_1, eta_ap, v_p_1, t_p_1, theta_ap_1,
               n_step=1000):
    # Initialize the alpha-proton charge and mass ratios.

    z_a = 2.
    mu_a = 4.

    # Initialise.

    d_r = (r_0 - r_1) / (1. * n_step)

    r = r_1
    n_p = n_p_1
    v_p = v_p_1
    t_p = t_p_1
    theta_ap = theta_ap_1

    theta_ap_min = 0.01
    theta_ap_max = 25.

    try:
        is_list_like = True
        temp = eta_ap[0]
    except:
        is_list_like = False

    # Loop.
    save_theta = []
    save_radius =[]
    for i in range(n_step):

        r = r_1 + ((i + 1) * d_r)

        n_p = n_p_1 * (r / r_1) ** -1.8
        v_p = v_p_1 * (r / r_1) ** -0.2
        t_p = t_p_1 * (r / r_1) ** -0.77

        arg_ = ((n_p ** 0.5 / t_p ** 1.5) * (z_a * (mu_a + 1) / (theta_ap + mu_a)) *
                (1 + (z_a ** 2 * eta_ap / theta_ap)) ** 0.5)

        if arg_ == 0:
            arg_ = math.exp(9)
        else:
            pass

        lambda_ap = 9 - np.log(arg_)

        d_theta_ap = ((-2.60e7) * ((n_p / (v_p * t_p ** 1.5))) * (
                mu_a ** 0.5 * z_a ** 2 / (eta_ap + 1) ** 2.5) *
                      ((theta_ap - 1.) * (eta_ap * theta_ap + 1.) ** 2.5 / (
                              theta_ap + mu_a) ** 1.5) * (lambda_ap) *
                      (d_r))

        theta_ap = theta_ap + d_theta_ap
        save_theta.append(theta_ap)
        save_radius.append(r)

        if (is_list_like):
            tk_i = np.where(theta_ap < theta_ap_min)
            theta_ap[tk_i] = theta_ap_min
        else:
            theta_ap = max([theta_ap, theta_ap_min])

        if (is_list_like):
            tk_i = np.where(theta_ap > theta_ap_max)
            theta_ap[tk_i] = theta_ap_max
        else:
            theta_ap = min([theta_ap, theta_ap_max])

    fn1 = 'theta' + str(theta_ap_1) + '.txt'
    fn2 = 'radius' + str(theta_ap_1) + '.txt'
    #np.savetxt(fn1, save_theta)
    #np.savetxt(fn2, save_radius)

    return theta_ap


def theta_loop(
        time,
        wind_radius,
        psp_radius,
        density_p,
        density_ap,
        speed,
        temp,
        theta,
):

    L = len(time)
    final_theta = np.zeros(L)
    for i in range(L):
        final_theta[i] = theta_ap_0(wind_radius[i], psp_radius[i], density_p[i],
                                    density_ap[i], speed[i], temp[i], theta[i])
        print('\r', f"{(i / L) * 100:.2f} %", end="")

    return final_theta


def radius_split(
        solar_data,
        spc_data,
        scalar_temps,
):
    dp_number = const.dp_number

    solar_sorted_data = {}
    spc_sorted_data = {}
    temp_sorted_data = {}

    radial_size = int(1 / (10 ** -dp_number))

    for i in range(radial_size):
        arg_ = round(i * (10 ** -dp_number), (dp_number + 1))
        solar_sorted_data[arg_] = {}
        spc_sorted_data[arg_] = {}
        temp_sorted_data[arg_] = {}
        for particle in solar_data.keys():
            solar_sorted_data[arg_][particle] = {}
            for key in solar_data[particle].keys():
                solar_sorted_data[arg_][particle][key] = []
        for spacecraft in spc_data.keys():
            spc_sorted_data[arg_][spacecraft] = {}
            for key in spc_data[spacecraft].keys():
                spc_sorted_data[arg_][spacecraft][key] = []
        for temp in scalar_temps.keys():
            temp_sorted_data[arg_][temp] = []

    for i in range(len(spc_data[psp]['RADIAL_DISTANCE_AU'])):
        arg_ = round(spc_data[psp]['RADIAL_DISTANCE_AU'][i], dp_number)
        for radius in solar_sorted_data.keys():
            if arg_ == radius:

                for particle in solar_sorted_data[radius]:
                    for key in solar_sorted_data[radius][particle].keys():
                        solar_sorted_data[radius][particle][key].append(
                            solar_data[particle][key][i])

                for spacecraft in spc_sorted_data[radius]:
                    for key in spc_sorted_data[radius][spacecraft].keys():
                        spc_sorted_data[radius][spacecraft][key].append(
                            spc_data[spacecraft][key][i])
                for temp in temp_sorted_data[radius]:
                    temp_sorted_data[radius][temp].append(scalar_temps[temp][i])
            else:
                pass
        print(f"{(i / len(spc_data[psp]['RADIAL_DISTANCE_AU'])) * 100:.2f} %", end="\r")

    return solar_sorted_data, spc_sorted_data, temp_sorted_data


def make_theta_vals(
        solar_data,
        spc_data,
        psp_scalar_temps,
        wind_radius_,

):
    time = solar_data[p]['time']
    density_p = solar_data[p]['np1']
    density_ap = psp_scalar_temps['dens_ap']
    temp = psp_scalar_temps['proton_1_k']
    speed = solar_data[p]['v_mag']
    theta = psp_scalar_temps['theta_ap']

    wind_radius = np.full(shape=len(spc_data['Wind_Orbit.csv'][t]),
                          fill_value=wind_radius_,
                          dtype=float)
    psp_radius = spc_data[psp]['RADIAL_DISTANCE_AU']

    final_theta = theta_loop(time, wind_radius, psp_radius, density_p, density_ap, speed,
                             temp, theta)

    return final_theta
