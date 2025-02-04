import os
import warnings

from constants import *

from modules.core.time import tictoc as stopwatch, convert as converter
from modules.core.constants import initialise_constants
from modules.core.variables import num_man as nm, string_man as sm
from modules.core.features import graph as graph
from modules.core.loadsave import file_dir as fdir

from modules.collisions.loadsave import loadsave as rw
from modules.collisions.features import lat_lon as lat_lon, \
    scalar_generate as sc_gen, scrub as scrub
from modules.collisions.model import theta_ap as theta_ap, errors as error

initialise_constants()
slash = fdir.slash()

t = 'time'
p = 'proton'
a = 'alpha'

particle_list = [p, a]
valid_enc = [4, 6, 7]


def load_generate(
        encounter,
        radius,

):
    parent_path = sm.slash_check(fdir.dir_parent())
    parent_path = cm.remove_end(parent_path)
    name = fdir.dir_name(parent_path)

    l = len(name) + 1

    for i in range(l):
        parent_path = cm.remove_end(parent_path)

    path = sm.jwos(parent_path, slash, 'data', slash, 'save', slash)

    if encounter == 0:
        enc = 'EA'
    else:
        enc = sm.jwos('E', str(encounter))

    fdir.dir_make(enc, path)

    path = sm.jwos(path, enc, slash)
    radius_dir = fdir.dir_make(str(radius), path)

    if not radius_dir:
        h = 1
        while h > 0:
            user_overide = input("Indication that these files already exist, "
                                 "would you like to overide them? (y/n)")
            if not isinstance(user_overide, str):
                raise TypeError(
                    "Error: User input is invalid."
                )
            else:
                var = cm.lower_all_letter(user_overide)
                if var == 'y':
                    h = 0
                elif var == 'n':
                    exit()
                else:
                    raise ValueError(
                        "Error: Only a Yes, No answer in the form of "
                        "y/n is required."
                    )

    return sm.jwos(path, str(radius), slash)


def encounter_generator(
        wind_rad=1,
        error_file=False
):
    stopwatch.start_time()

    if not isinstance(wind_rad, (int, float)):
        raise TypeError(
            "Error: Wind radius must be either a float or integer, "
            f"instead got type of {type(wind_rad)}"
        )

    if not isinstance(error_file, bool):
        raise TypeError(
            "Error:"
        )

    print('Currently loaded encounters:', valid_enc, '\n')

    valid_full = ['f', 'F', 'full', 'Full', 'ful', 'Ful', 'ff', 'FF', 'fF', 'Ff']
    valid_single = ['s', 'S', 'single', 'Single', 'ss', 'SS', 'sS', 'Ss']

    h = 1
    while h > 0:
        data_set_input = input('Full data set or singular? (F/S)')
        if data_set_input in valid_full:
            encount = 0
            h = 0
        elif data_set_input in valid_single:
            h = 0
            g = 1
            while g > 0:
                enc_input = input('Please enter an encounter:')
                if sm.valid_string(enc_input):
                    if int(enc_input) in valid_enc:
                        encount = int(enc_input)
                        g = 0
                    elif enc_input == '':
                        print('Error: No input provided.')
                    else:
                        print('Error: No corresponding encounter available.')
                else:
                    print('Error: Argument provided is not valid,'
                          f' argument {enc_input}, is of type {type(enc_input)}.')
        else:
            print('Error: Please make a valid selection.')

    # Generate directory strings for encounters
    encounter_number = encount
    enc(encount, valid_enc, error_file)
    print('\n')
    error_files = enc.error_files_loaded
    mm_data = rw.encounter_import(encount, valid_enc, error_files)
    if error_files:
        error_data = rw.error_import()
    sc_data = rw.sc_import()
    save_loc = load_generate(encount, wind_rad)
    print('\nData import sucessful.\n')

    # Generate equal sized arays for all data

    for x in enc.encounter:
        res = list(mm_data[x].keys())[0]
        t_ = mm_data[x][res][t]

        for y in mm_data[x].keys():
            xp = mm_data[x][y][t]
            for z in mm_data[x][y].keys():
                fp = mm_data[x][y][z]
                mm_data[x][y][z] = np.interp(t_, xp, fp)
                # mm_data[x][y][z] = np.resize(mm_data[x][y][z], min_len)

        if error_files:
            for y in error_data[x].keys():
                xp = error_data[x][y][t]
                for z in error_data[x][y].keys():
                    fp = error_data[x][y][z]
                    error_data[x][y][z] = np.interp(t_, xp, fp)
                    # error_data[x][y][z] = np.resize(error_data[x][y][z], min_len)

        t_
        for y in sc_data[x].keys():
            xp = sc_data[x][y][t]
            for z in sc_data[x][y].keys():
                fp = sc_data[x][y][z]
                sc_data[x][y][z] = np.interp(t_, xp, fp)
                # sc_data[x][y][z] = np.resize(sc_data[x][y][z], min_len)

    # h5g files

    # Generate temperatures and velocity magnitudes
    print('Generating velocity magnitudes and temperature file... \n')
    mm_data = sc_gen.scalar_velocity(mm_data)

    psp_temps, wind_temps = sc_gen.scalar_temps(mm_data, sc_data)


    #tol_value = 5
    #guess = {4: 11.5, 6: 3.5, 7: 7.8}
    #if encounter_number == 0:
    #    for key in psp_temps.keys():
    #        value = int(cm.remove_begin(key))
    #        print(key, value)
    #        print(type(psp_temps[key]['theta_ap']))
    #        psp_temps[key]['theta_ap'] = theta_ap.remove_theta(psp_temps[key]['theta_ap'],
    #                                                           guess[value], tol_value)
    #else:
    #    arg_ = str('E' + str(encounter_number))
    #    psp_temps[arg_]['theta_ap'] = theta_ap.remove_theta(psp_temps[arg_]['theta_ap'],
    #                                                        guess[encounter_number],
    #                                                        tol_value)

    # Generate file and/or combine files (remember to do the scalar temp files)
    print('Generating data file... \n')
    solar_data = {}
    errors = {}
    spc_data = {}
    psp_scalar_temps = {}
    wind_scalar_temps = {}

    for particle in particle_list:
        solar_data[particle] = {}
        errors[particle] = {}

    for key in psp_temps.keys():
        for value in psp_temps[key].keys():
            psp_scalar_temps[value] = []

    for key in wind_temps.keys():
        for value in wind_temps[key].keys():
            wind_scalar_temps[value] = []

    for key in psp_scalar_temps.keys():
        for encount in psp_temps.keys():
            for i in range(len(psp_temps[encount][key])):
                psp_scalar_temps[key].append(psp_temps[encount][key][i])

    for key in wind_scalar_temps.keys():
        for encount in wind_temps.keys():
            for i in range(len(wind_temps[encount][key])):
                wind_scalar_temps[key].append(wind_temps[encount][key][i])

    for x in range(1):
        encount = enc.encounter[x]
        for y in range(1):
            for particle in particle_list:
                indx = particle_list.index(particle)
                if nm.is_even(indx):
                    arg_x_ = 0
                elif nm.is_odd(indx):
                    arg_x_ = 1
                for z in mm_data[encount][enc.encounter_names[2 * x + arg_x_]].keys():
                    solar_data[particle][z] = []

        if error_files:
            for y in range(1):
                for particle in particle_list:
                    indx = particle_list.index(particle)
                    if nm.is_even((indx)):
                        arg_x_ = 0
                    else:
                        arg_x_ = 1
                    for z in error_data[encount][
                        enc.encounter_errors[2 * x + arg_x_]].keys():
                        errors[particle][z] = []

        for y in range(len(enc.sc_names)):
            spc_data[enc.sc_names[y]] = {}
            for z in sc_data[encount][enc.sc_names[y]].keys():
                spc_data[enc.sc_names[y]][z] = []

    for x in range(enc.num_of_encs):
        encount = enc.encounter[x]
        for particle in particle_list:
            indx = particle_list.index(particle)
            if nm.is_even(indx):
                arg_x_ = 0
            elif nm.is_odd(indx):
                arg_x_ = 1
            for z in solar_data[particle].keys():
                for w in range(
                        len(mm_data[encount][enc.encounter_names[2 * x + arg_x_]][z])):
                    solar_data[particle][z].append(
                        mm_data[encount][enc.encounter_names[2 * x + arg_x_]][z][w])
            if error_files:
                for z in errors[particle].keys():
                    for w in range(
                            len(error_data[encount][enc.encounter_errors[2 * x + arg_x_]][
                                    z])):
                        errors[particle][z].append(
                            error_data[encount][enc.encounter_errors[2 * x + arg_x_]][z][
                                w])

        for y in enc.sc_names:
            for z in spc_data[y].keys():
                for w in range(len(sc_data[encount][y][z])):
                    spc_data[y][z].append(sc_data[encount][y][z][w])

    print('Scrubbing data...')
    solar_data, errors, spc_data = scrub.scrub_data(solar_data, errors, spc_data, error_files)

    spc_data[enc.sc_names[2]] = lat_lon.latlong_psp(spc_data[enc.sc_names[2]])
    spc_data[enc.sc_names[1]] = lat_lon.latlong_wind(spc_data[enc.sc_names[1]])

    # Generate single time set for the whole data set in appropriate unit
    solar_data[t] = []
    for i in range(len(solar_data[p][t])):
        solar_data[t].append(converter.epoch_time(solar_data[p][t][i]))


    theta_ap_0 = psp_scalar_temps['theta_ap']
    wind_radius = wind_rad
    theta_ap_final = theta_ap.make_theta_vals(solar_data, spc_data, psp_scalar_temps,
                                              wind_radius)
    print('Note: Files have been generated and loaded in.', '\n')

    if error_files:
        uncertain = error.sigma_value(solar_data, errors, psp_scalar_temps)

    theta = {'0.1 - 0.2': theta_ap_0, str(wind_rad): theta_ap_final}

    print(save_loc)
    file_names = ['theta_i.txt', 'theta_f.txt', 'wind_theta.txt']
    temp_files = [theta['0.1 - 0.2'], theta[str(wind_rad)],
                  wind_scalar_temps['wind_theta']]

    if len(file_names) != len(temp_files):
        warnings.warn("Warning: File save system has unequal lengths.")
    else:
        for i in range(len(file_names)):
            loc = sm.jwos(save_loc, file_names[i])
            np.savetxt(loc, temp_files[i])

    return



def uncertain(

):
    x = np.linspace(0, 15, 1000)
    y = uncertain[p]['Isotropy']
    z = uncertain[a]['Isotropy']
    print(y, len(uncertain[p]['Scalar Temp']))
    print(z)

    import math
    from scipy.optimize import curve_fit
    def maxwellian(x, r, m, s):
        return (r / (s * np.sqrt(math.pi))) * np.exp(- (x - m) ** 2 / (2 * (s ** 2)))

    y = y[np.logical_not(np.isnan(y))]
    z = z[np.logical_not(np.isnan(z))]
    bn_i = 10
    bn_a = 5
    data_entries, bins = np.histogram(y, bins=bn_i, density=True)
    data_alpha, bins_alpha = np.histogram(z, bins=bn_a, density=True)
    binscenters = np.array([0.5 * (bins[i] + bins[i + 1]) for i in range(len(bins) - 1)])
    binalphacenters = np.array(
        [0.5 * (bins_alpha[i] + bins_alpha[i + 1]) for i in range(len(bins_alpha) - 1)])
    from modules.core.features import smooth as smoothing
    popt, pcov = curve_fit(maxwellian, xdata=binscenters, ydata=data_entries,
                           maxfev=10000)
    popta, pcova = curve_fit(maxwellian, xdata=binalphacenters, ydata=data_alpha,
                             maxfev=10000)
    xspace = np.linspace(0, 20, 10000)
    yspace = maxwellian(xspace, *popt)
    x = np.linspace(0, 20, 10000)
    zspace = maxwellian(x, *popta)

    data = {'Proton': yspace, 'Alpha': zspace}
    print(zspace)
    color = ['black', 'blue']
    graph.histogram(xspace, y)
    graph.histogram(xspace, z, )
    graph.graph(xspace, data, colours=color, title='',
                x_axis=r'$\frac{\sigma_{|v|}}{|v|}$', y_axis='Probability Density',
                limits=False, x_lim=25, y_lim=5)

    intervals_per_day = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 20, 23, 48,
                         96, 192]
    interval_length = [24. / i for i in intervals_per_day]
    vals_ = {}
    vals_[p] = []
    vals_[a] = []

    for i in range(len(intervals_per_day)):
        e_p_val, e_a_val = error.gen_uncer(solar_data, psp_scalar_temps,
                                           intervals_per_day[i])
        vals_[p].append(e_p_val)
        vals_[a].append(e_a_val)

    graph(interval_length, vals_, title='', x_axis='Interval length hours',
          y_axis='Median Std', log=False)  # interval length log

    stopwatch.end_time()

    return solar_data, errors, psp_scalar_temps, wind_scalar_temps, vals_
