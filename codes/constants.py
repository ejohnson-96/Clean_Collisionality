import numpy as np
from core.rw import file_dir as fd


def file_dir_gen(
        enc,
        valid_enc,
):
    path = fd.dir_path()
    str_load = str(path + "/data/load/")
    str_save = str(path + "/data/save/")

    if enc == 0:
        L = len(valid_enc)
        encounter = []
        for i in range(L):
            encounter[i] = valid_enc[i]
            encounter.append('E' + str(int(encounter[i])))
    else:
        L = 1
        encounter = []
        encounter.append('E' + str(enc))
    num_of_encs = L
    encounter_names = []
    encounter_errors = []

    for i in range(L):
        val = (encounter[i])
        encounter_names.append(val + '_protons.csv')
        encounter_names.append(val + '_alphas.csv')
        encounter_errors.append(val + '_proton_errors.csv')
        encounter_errors.append(val + '_alpha_errors.csv')

    sc_names = []
    sc_names.append('PSP.csv')
    sc_names.append('Wind_Orbit.csv')
    sc_names.append('PSP_Orbit.csv')
    sc_names.append('Wind_Outside_Range_Hour.csv')
    sc_names.append('Wind_Outside_Range_Min.csv')
    sc_names.append('Wind_Temps.csv')

    num_of_sc = len(sc_names)


    return


def constants(

):


    return