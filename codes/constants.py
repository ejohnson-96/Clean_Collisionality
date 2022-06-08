import numpy as np
from core.rw import file_dir as fd


def enc(
        encount,
        valid_enc,
        error_files=False,
):
    enc.path = fd.dir_path()
    enc.str_load = str(enc.path + "/data/load/")
    enc.str_save = str(enc.path + "/data/save/")
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
            val_ = (enc.encounter[i])
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

    return


def constants(

):


    return


def initialise_constants(

):
    constants()

    return