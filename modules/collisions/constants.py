import numpy as np
from modules.core.loadsave import file_dir as fd
from modules.core.system import config as sys_con
from modules.core.variables import char_man as cm

slash = fd.slash()

def enc(
        encount,
        valid_enc,
        error_files=False,
):
    path = fd.dir_parent()
    for i in range(7):
        path = cm.remove_end(path)
    enc.path = path
    enc.str_load = str(enc.path + "data" + slash + "load" + slash)
    enc.str_save = str(enc.path + "data" + slash + "save" + slash)
    enc.error_files_loaded = error_files

    if encount == 0:
        L = len(valid_enc)
        enc.encounter = []
        for i in range(L):
            arg_ = valid_enc[i]
            enc.encounter.append('E' + str(int(arg_)))
    else:
        L = 1
        enc.encounter = []
        enc.encounter.append('E' + str(encount))
    enc.num_of_encs = L
    enc.encounter_names = []
    enc.encounter_errors = []

    for i in range(L):
        val_ = (enc.encounter[i])
        enc.encounter_names.append(val_ + '_protons.csv')
        enc.encounter_names.append(val_ + '_alphas.csv')
    if error_files:
        for i in range(L):
            enc.encounter_errors.append(val_ + '_proton_errors.csv')
            enc.encounter_errors.append(val_ + '_alpha_errors.csv')

    enc.sc_names = []
    enc.sc_names.append('PSP.csv')
    enc.sc_names.append('Wind_Orbit.csv')
    enc.sc_names.append('PSP_Orbit.csv')
    enc.sc_names.append('Wind_Outside_Range_Hour.csv')
    enc.sc_names.append('Wind_Outside_Range_Min.csv')
    enc.sc_names.append('Wind_Temps.csv')

    enc.num_of_sc = len(enc.sc_names)
    if error_files:
        arg_errors_ = len(enc.encounter_errors)
    else:
        arg_errors_ = 0
    enc.num_files = len(enc.encounter_names) + arg_errors_ + len(enc.sc_names)

    return

def const(

):
    # Min and Max Values#
    density_max = 10 ** 25
    density_min = 10 ** -6

    speed_max = 10 ** 3
    speed_min = 0

    b_field_max = 10 ** 8
    b_field_min = 0

    temp_max = 10 ** 6
    temp_min = 0

    chi_squared_max = 5000
    chi_squared_min = 0

    dr_max = 1
    dr_min = 0

    dv_max = 10
    dv_min = 0

    dT_perp_max = 10
    dT_perp_min = 0

    bmag_max = 9900
    bmag_min = 0

    v_gse_max = 99900
    v_gse_min = 0

    dens_add_max = 990
    dens_add_min = 0

    nanp_max = 9
    nanp_min = 0

    presure_max = 90
    presure_min = 0

    c_max = 10 ** 6
    c_min = 0

    const.data_units = {}

    const.data_units[0] = [0, 1, 1, 2, 2, 2, 3, 3, 3, 2, 4, 4, 4, 4, 2]
    const.data_units[1] = [0, 1, 2, 2, 2, 3, 3, 3, 4, 4, 2]

    const.error_units = {}

    const.error_units[0] = [0, 5, 1, 8, 8, 5, 5, 5, 1, 8, 6, 7]
    const.error_units[1] = [0, 5, 1, 8, 8, 5, 5, 5, 1, 8, 6, 7]

    const.sc_units = {}

    const.sc_units[0] = [0, 0, 0, 0, 3, 3, 3, 3, 2, 2, 2, 2, 0, 0, 0, 0]
    const.sc_units[1] = [0, 0, 0, 0, 3, 3, 3, 3, 2, 0, 0, 0, 0]

    const.sc_units[2] = [0, 0, 0, 0, 0, 0]
    const.sc_units[3] = const.sc_units[1]
    const.sc_units[4] = [0, 9, 9, 9, 9, 9, 9, 9, 9, 2, 10, 10, 2, 11, 14, 12, 13, 11, 11,
                         11]
    const.sc_units[5] = [0, 1, 4, 1, 4]

    const.var_max = [(10 ** 30), density_max, speed_max, b_field_max, temp_max,
                     chi_squared_max, dr_max, dv_max, dT_perp_max, bmag_max, v_gse_max,
                     dens_add_max, nanp_max, presure_max, c_max]

    const.var_min = [0, density_min, speed_min, b_field_min, temp_min, chi_squared_min,
                     dr_min, dv_min, dT_perp_min, bmag_min, v_gse_min, dens_add_min,
                     nanp_min, presure_min, c_min]


    return

def load_constants(

):
    return
