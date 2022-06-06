import os
import pathlib


def types(

):
    types.valid_file_types = ['csv', 'xlsx', 'pkl']
    types.valid_read_types = ['_csv', '_excel', '_pickle']

    res = 'All valid file types loaded.'

    return res


def dir_path(

):
    directory = os.getcwd()
    path = str(pathlib.Path(directory).parent)

    return path
