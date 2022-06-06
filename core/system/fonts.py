import os
from matplotlib import font_manager
from core.variables import char_man as cm
from core.variables import string_man as sm


def list_fonts(

):
    fpaths = font_manager.findSystemFonts()

    res = {'name': [], 'dir': []}

    i = 0
    for font_name in fpaths:
        f = font_manager.get_font(font_name)
        res['name'].append(f.family_name)
        res['dir'].append(str(fpaths[i]))
        i += 1

    return res


def font_names(

):
    fonts = list_fonts()
    res = fonts['name']

    return res


def font_dir(

):
    fonts = list_fonts()
    res = fonts['dir']

    return res


def check_font_valid(
        font,
):
    if not isinstance(font, str):
        raise ValueError(
            'Error: Argument must be a string, instead '
            f'got type {type(font)}'
        )

    valid_fonts = font_names()

    for i in range(len(valid_fonts)):
        valid_fonts[i] = cm.lower_all_letter(valid_fonts[i])

    input_font = cm.lower_all_letter(font)

    if input_font in valid_fonts:
        res = True
    else:
        res = False

    return res


def load_font(
        loc,
):
    if not isinstance(loc, str):
        raise TypeError(
            'Error: Argument must be a string, instead '
            f'got type {type(loc)}'
        )

    valid_file = 'ttf'

    file_check = loc.split('.')

    ext = cm.lower_all_letter(file_check[-1])

    if ext == valid_file:
        font_manager.fontManager.addfont(loc)
    else:
        raise ValueError(
            'Error: Directory location provided does not appear'
            ' to be the correct file type.'
        )

    return
