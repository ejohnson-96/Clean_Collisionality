import random


def hex_to_rgb(
        hex,
):
    if not isinstance(hex, str):
        raise TypeError(
            'Error: Argument must be a string, instead'
            f' got type {type(hex)}.'
        )

    if '#' not in hex:
        raise ValueError(
            'Error: String must begin with #.'
        )
    else:
        h = hex.lstrip('#')

    return tuple(int(h[i:i+2], 16) for i in (0, 2, 4))

def rgb_to_hex(
        r,
        g,
        b,
):
    for val in (r,g,b):
        if not isinstance(val, int):
            raise TypeError(
                'Error: '
            )
    rgb = (r,g,b)

    return '#%02x%02x%02x' % rgb


def rand_rgb(

):
    r = random.random()
    g = random.random()
    b = random.random()

    return r, g, b