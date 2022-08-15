import sys
from modules.core.variables.string_man import valid_string

windows = 'Windows'
mac = 'Mac OS'
linux = 'Linux'


def get_platform():
    platforms = {
        'linux1': linux,
        'linux2': linux,
        'darwin': mac,
        'win32': windows
    }
    if sys.platform not in platforms:
        return sys.platform

    return platforms[sys.platform]


def windows_os(

):
    if get_platform() == windows:
        return True
    else:
        return False


def mac_os(

):
    if get_platform() == mac:
        return True
    else:
        return False


def linux_os(

):
    if get_platform() == linux:
        return True
    else:
        return False

