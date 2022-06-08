import numpy as np
from modules.core.rw import file_dir as fd
from modules.core.system import config as sys_con
from modules.core.variables import char_man as cm

windows_check = sys_con.windows_os()
if windows_check:
    slash = '\\'
else:
    slash = '/'


def enc(
        encount,
        valid_enc,
        error_files=False,
):
    path = fd.dir_path()
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
            enc.encounter[i] = valid_enc[i]
            enc.encounter.append('E' + str(int(enc.encounter[i])))
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

    num_of_sc = len(enc.sc_names)

    return


def constants(

):
    return
