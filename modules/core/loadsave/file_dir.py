import os
import pathlib
import warnings
from modules.core.variables import string_man as sm
from modules.core.system import config as sys_con


def dir_parent(

):
    directory = os.getcwd()
    path = str(pathlib.Path(directory).parent)

    return path


def dir_path(

):
    directory = os.getcwd()
    path = str(pathlib.Path(directory))

    return path


path_dir = dir_path()


def dir_make(
        name,
        loc=path_dir,
):
    if not isinstance(name, str):
        raise TypeError(
            f"Directory name passed is not a string, "
            f"instead got {type(name)}"
        )

    path = sm.slash_check(loc) + name
    isExist = os.path.exists(path)

    if isExist:
        warnings.warn("Warning: Directory already exits.")
        return False
    else:
        os.mkdir(path)
        return True


def dir_name(
        loc=None,
):
    if isinstance(loc, type(None)):
        path = os.getcwd()
        return os.path.basename(path)
    else:
        if not isinstance(loc, str):
            raise TypeError(
                "Error: Directory path passed must be a string,"
                f" instead got type {type(loc)}."
            )
        isExist = os.path.exists(loc)
        if isExist:
            return os.path.basename(loc)
        else:
            warnings.warn(
                "Error: Directory location provided does not exist,"
                f" location provided is: {loc}"
            )
            return False


def file_list(
        loc=path_dir,
):
    path = sm.slash_check(loc)
    res = next(os.walk(path))[2]

    return res


def file_num(
        loc=path_dir,
):
    L = file_list(loc)
    res = len(L)

    return res


def folder_list(
        loc=path_dir,
):
    path = sm.slash_check(loc)
    res = next(os.walk(path))[1]

    return res


def folder_num(
        loc=path_dir,
):
    L = folder_list(loc)
    res = len(L)

    return res


def dir_list(
        loc=path_dir,
):
    path = sm.slash_check(loc)
    res = os.listdir(path)

    return res


def dir_num(
        loc=path_dir,
):
    L = dir_list(loc)
    res = len(L)

    return res


def slash(

):
    windows_check = sys_con.windows_os()

    if windows_check:
        return '\\'
    else:
        return '/'
