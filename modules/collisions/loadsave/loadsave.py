from modules.collisions.constants import *
from modules.core.rw import file_import as fi, file_dir as fd
from modules.core.system import config as sys_con

slash = fd.slash()

def create_encs(
        encount,
        valid_enc,
        error_files,
):
    if not isinstance(encount, int):
        raise TypeError(
            "Error: The encounter argument must be an integer,"
            f"instead got type {type(encount)}."
        )
    if not isinstance(error_files, bool):
        raise TypeError(
            "Error: Indicate if the error files are present "
            f"with a boolean argument, instead got {type(error_files)}."
        )
    if not isinstance(valid_enc, list):
        raise TypeError(
            "Error: The valid encounter must be a list "
            f" instead got {type(encount)}."
        )


    enc(encount, valid_enc, error_files)

    return



def encounter_import(
        encount,
        valid_enc,
        error_files,
):
    create_encs(encount, valid_enc, error_files)

    files = {}

    for i in range(enc.num_of_encs):
        print('Data File: ' + enc.encounter[i])
        files[enc.encounter[i]] = {}
        for j in range(2):
            val = str(enc.encounter[i] + slash + enc.encounter_names[j + 2 * i])
            files[enc.encounter[i]][enc.encounter_names[j + 2 * i]] = fi.file_import(
                enc.str_load + val)
    return files


def sc_import(

):
    files = {}

    for i in range(enc.num_of_encs):
        print('Spacecraft File: ' + enc.encounter[i])
        files[enc.encounter[i]] = {}
        for key in enc.sc_names:
            val = str(enc.encounter[i] + slash + 'Position' + slash + key)
            files[enc.encounter[i]][key] = fi.file_import(enc.str_load + val)

    return files



def error_import(

):
    files = {}

    if enc.error_files_loaded:
        for i in range(enc.num_of_encs):
            print('Error File: ' + enc.encounter[i])
            files[enc.encounter[i]] = {}
            for j in range(2):
                val = str(enc.encounter[i] + slash + enc.encounter_errors[j + 2 * i])
                files[enc.encounter[i]][enc.encounter_errors[j + 2 * i]] = fi.file_import(enc.str_load + val)
    else:
        raise ValueError(
            'Error: Tried to load error files when constant indicates none are present.'
        )

    return files
