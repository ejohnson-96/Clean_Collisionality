from core.rw import file_dir as fd


def file_dir_gen(
        enc,
        valid_enc,
):
    path = fd.dir_path()
    print(path)
    const.str_dir = str(path + "/data/load/")
    const.str_save = str(path + "/data/save/")

    if enc == 0:
        L = len(valid_enc)
        encounter = np.zeros(L)
        const.encounter = []
        for i in range(L):
            encounter[i] = valid_enc[i]
            const.encounter.append('E' + str(int(encounter[i])))
    else:
        L = 1
        encounter = np.zeros(L)
        encounter[0] = enc
        const.encounter = []
        const.encounter.append('E' + str(int(encounter[0])))
    const.num_of_encs = L
    const.encounter_names = []
    const.encounter_errors = []

    for i in range(L):
        val = (const.encounter[i])
        const.encounter_names.append(val + '_protons.csv')
        const.encounter_names.append(val + '_alphas.csv')
        const.encounter_errors.append(val + '_proton_errors.csv')
        const.encounter_errors.append(val + '_alpha_errors.csv')

    const.sc_names = []
    const.sc_names.append('PSP.csv')
    const.sc_names.append('Wind_Orbit.csv')
    const.sc_names.append('PSP_Orbit.csv')
    const.sc_names.append('Wind_Outside_Range_Hour.csv')
    const.sc_names.append('Wind_Outside_Range_Min.csv')
    const.sc_names.append('Wind_Temps.csv')

    const.num_of_sc = len(const.sc_names)


    return


def constants(

):

    file_dir_gen()

    return