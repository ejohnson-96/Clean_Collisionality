import keyboard
from modules.core.variables import char_man as cm


def user_input(

):

    key = keyboard.is_pressed()

    return str(cm.lower_all_letter(key))


def navigation(
        key_press,
        north=0,
        east=0,
        south=0,
        west=0,
):

    if not isinstance(key_press, str):
        raise ValueError(
            f'Error: Key not found, got {key_press}.'
        )

    for val in (north, east, south, west):
        if not isinstance(val, int):
            raise ValueError(
                'Error: Co-ordinates must be integers, '
                f'instead got {type(val)}.'
            )

    if key_press == 'w':
        north = north + 1
    if key_press == 'd':
        east = east + 1
    if key_press == 's':
        south = south + 1
    if key_press == 'a':
        west = west + 1

    return key_press

