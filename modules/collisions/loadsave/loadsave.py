from modules.collisions.constants import *
from modules.core.rw import file_import as fi


def create_encs(
        encount,
        valid_enc,
):
    if not isinstance(encount, int):
        raise TypeError(
            "Error: The encounter argument must be an integer,"
            f"instead got type {type(encount)}."
        )

    enc(encount, valid_enc)

    return



def encounter_import(
        encount,
        valid_enc,

):
    create_encs(encount, valid_enc)

    files = {}

    for i in range(enc.num_of_encs):
        print('Data File: ' + enc.encounter[i])
        files[enc.encounter[i]] = {}
        for j in range(2):
            val = str(enc.encounter[i] + '/' + enc.encounter_names[j + 2 * i])
            files[enc.encounter[i]][enc.encounter_names[j + 2 * i]] = fi.file_import(
                enc.str_load + val)
    return files


def sc_import(

):
    files = {}

    for i in range(enc.num_of_encs):
        print('Spacecraft for ' + enc.encounter[i])
        files[enc.encounter[i]] = {}
        for key in enc.sc_names:
            val = str(enc.encounter[i] + '/Position/' + key)
            files[enc.encounter[i]][key] = fi.file_import(enc.str_load + val)
    print('\n',
          'Warning: Please ensure all data is in the correct time range for the encounter.',
          '\n')
    return files



def error_import(

):
    files = {}

    for i in range(const.num_of_encs):
        print('Error File: ' + const.encounter[i])
        files[const.encounter[i]] = {}
        for j in range(2):
            val = str(const.encounter[i] + '/' + const.encounter_errors[j + 2 * i])
            files[const.encounter[i]][
                const.encounter_errors[j + 2 * i]] = fi.file_import(const.str_load + val)
    return files
