import sys

windows = 'Windows'
mac  = 'Mac OS'
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
    os = get_platform()

    if os == windows:
        return True
    else:
        return False


def mac_os(

):
    os = get_platform()

    if os == mac:
        return True
    else:
        return False


def linux_os(

):
    os = get_platform()

    if os == linux:
        return True
    else:
        return False