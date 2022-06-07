from modules.collisions.constants import enc as enc
from modules.core.rw import file_import as fi


def encounter_import(

):
    files = {}

    for i in range(enc.num_of_encs):
        print('Data File: ' + enc.encounter[i])
        files[enc.encounter[i]] = {}
        for j in range(2):
            val = str(enc.encounter[i] + '/' + enc.encounter_names[j + 2 * i])
            files[enc.encounter[i]][enc.encounter_names[j + 2 * i]] = fi.file_import(
                enc.str_dir + val)
    return files


def sc_import(

):
    files = {}

    for i in range(enc.num_of_encs):
        print('Spacecraft for ' + enc.encounter[i])
        files[enc.encounter[i]] = {}
        for key in enc.sc_names:
            val = str(enc.encounter[i] + '/Position/' + key)
            files[enc.encounter[i]][key] = fi.file_import(enc.str_dir + val)
    print('\n',
          'Warning: Please ensure all data is in the correct time range for the encounter.',
          '\n')
    return files


def epoch_time(
        epoch_time,
):
    import datetime
    res = datetime.datetime.fromtimestamp(epoch_time)

    return res


def error_import(

):
    files = {}

    for i in range(const.num_of_encs):
        print('Error File: ' + const.encounter[i])
        files[const.encounter[i]] = {}
        for j in range(2):
            val = str(const.encounter[i] + '/' + const.encounter_errors[j + 2 * i])
            files[const.encounter[i]][
                const.encounter_errors[j + 2 * i]] = fimp.file_import(const.str_dir + val)
    return files
