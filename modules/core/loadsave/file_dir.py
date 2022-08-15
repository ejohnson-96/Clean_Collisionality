import os
import pathlib
import warnings
from modules.core.variables import string_man as sm
from modules.core.system import config as sys_con


def dir_parent(

):
    return str(pathlib.Path(os.getcwd()).parent)


def dir_path(

):
    return str(pathlib.Path(os.getcwd()))


path_dir = dir_path()


def dir_make(
        name,
        loc=path_dir,
):
    if not isinstance(name, str):
        raise TypeError(
            f"Directory name {name} passed is not a string, "
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
        return os.path.basename(os.getcwd())
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
    return next(os.walk(sm.slash_check(loc)))[2]


def file_num(
        loc=path_dir,
):
    return len(file_list(loc))


def folder_list(
        loc=path_dir,
):
    return next(os.walk(sm.slash_check(loc)))[1]


def folder_num(
        loc=path_dir,
):
    return len(folder_list(loc))


def dir_list(
        loc=path_dir,
):
    return os.listdir(sm.slash_check(loc))


def dir_num(
        loc=path_dir,
):
    return len(dir_list(loc))


def slash(

):
    if sys_con.windows_os():
        return '\\'
    else:
        return '/'


def file_ext(
        loc=path_dir,
):
    if not isinstance(loc, str):
        raise TypeError(
            'Error: File location must be of type string,'
            f' instead got type of {type(loc)}.'
        )

    return os.path.splitext(loc)[-1].lower()


def all_file_type(
        file_type,
        loc=path_dir,
):
    loc = sm.slash_check(loc)
    res = []
    for root, dirs, files in os.walk(loc):
        for file in files:
            if file.endswith(file_type):
                res.append(os.path.join(root, file))
    return res
