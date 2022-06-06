import os
import ctypes
import random
import string

val_sym = "!@#$%^&*()"

def is_user_admin(

):
    try:
        res = (os.getuid() == 0)
    except AttributeError:
        res = ctypes.windll.shell32.IsUserAnAdmin() != 0
    return res


def list_users(

):

    return

def pwd_gen(
        len=16,
        numbers=True,
        symbols=True,
):
    min = 8

    if not isinstance(len, int):
        raise TypeError(
            "Error: Argument must be an integer."
        )

    if not numbers:
        x = ''
    else:
        x = string.digits

    if not symbols:
        y = ''
    else:
        y = val_sym

    characters = list(string.ascii_letters + x + y)




    if len < min:
        raise ValueError(
            f"Error: Integer must be greater than {min} characters."
        )
    else:
        password = []
        for i in range(len):
            password.append(random.choice(characters))
        random.shuffle(password)

        res = "".join(password)

    return res


#import wmi
def os_vers(

):
    #try:

    #except:
     #   data = wmi.WMI()
      #  for os_name in data.Win32_OperatingSystem():
       #     print(os_name.Caption)

    return