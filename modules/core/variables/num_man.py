import math
import random


def valid_num(
        num,
):
    if not isinstance(num, (int, float)):
        raise TypeError(
            f"Error: Argument {num} is not of type int or float, "
            f"instead got type of {type(num)}."
        )
    return num


def valid_int(
        num
):
    if not isinstance(num, int):
        raise TypeError(
            f"Error: Argument {num} is not of type int or float, "
            f"instead got type of {type(num)}."
        )
    else:
        return num


def valid_float(
        num
):
    if not isinstance(num, (int, float)):
        raise TypeError(
            f"Error: Argument {num} is not of type int or float, "
            f"instead got type of {type(num)}."
        )
    return num


def add(
        a,
        b,
        c=0,
        d=0,
        e=0,
        f=0,
        g=0,
        h=0,
        i=0,
        j=0,
        k=0,
        l=0,
        m=0,
):
    for arg_name in (a, b, c, d, e, f, g, h, i, j, k, l, m):
        arg_name = valid_num(arg_name)

    return a + b + c + d + e + f + g + h + i + j + k + l + m


def sub(
        a,
        b,
        c=0,
        d=0,
        e=0,
        f=0,
        g=0,
        h=0,
        i=0,
        j=0,
        k=0,
        l=0,
        m=0,
):
    for arg_name in (a, b, c, d, e, f, g, h, i, j, k, l, m):
        arg_name = valid_num(arg_name)

    return a - b - c - d - e - f - g - h - i - j - k - l - m


def multi(
        a,
        b,
        c=0,
        d=0,
        e=0,
        f=0,
        g=0,
        h=0,
        i=0,
        j=0,
        k=0,
        l=0,
        m=0,
):
    for arg_name in (a, b, c, d, e, f, g, h, i, j, k, l, m):
        if not isinstance(arg_name, (int, float)):
            raise TypeError(
                f"Argument '{arg_name}'must be an integer or float."
            )

    arg_ = a * b

    for default_val in (c, d, e, f, g, h, i, j, k, l, m):
        if default_val == 0:
            break
        else:
            arg_ = arg_ * default_val

    res = arg_

    return res


def div(
        a,
        b,
):
    a = valid_num(a)
    b = valid_num(b)

    if b == 0:
        raise ValueError(
            "Error: Second argument cannot be zero dummy."
        )

    return a / b


def dot_product(
        a,
        b,
        theta=False,
        radians=False
):
    a = valid_num(a)
    b = valid_num(b)

    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if not isinstance(theta, bool):
            if not isinstance(theta, (int, float)):
                raise TypeError(
                    "Theta value must be an integer or float, instead"
                    f" got type {type(theta)}."
                )
        else:
            raise TypeError(
                "Theta value was not provided."
            )
        if radians:
            res = a * b * math.cos(theta)
        else:
            res = a * b * math.cos(deg_rad(theta))
    elif isinstance(a, (tuple, list)) and isinstance(b, (tuple, list)):
        arg_ = int((len(a) + len(b)) / 2)
        if arg_ == 3:
            res = a[0] * a[1] * a[2] + b[0] * b[1] * b[2]
        elif arg_ == 2:
            res = a[0] * a[1] + b[0] * b[1]
        else:
            raise ValueError(
                "Vector must have a dimension of either 2 or 3."
            )

    return res


def cross_product(
        a,
        b,
        theta=False,
        radians=False
):
    for arg_name in (a, b):
        if not isinstance(arg_name, (int, float, tuple, list)):
            raise TypeError(
                f"Argument '{arg_name}'must be an integer, float or"
                f" a vector."
            )
    if isinstance(a, (int, float)) and isinstance(b, (int, float)):
        if not isinstance(theta, bool):
            if not isinstance(theta, (int, float)):
                raise TypeError(
                    "Theta value must be an integer or float, instead"
                    f" got type {type(theta)}."
                )
        else:
            raise TypeError(
                "Theta value was not provided."
            )
        if radians:
            res = a * b * math.sin(theta)
        else:
            res = a * b * math.sin(deg_rad(theta))
    elif isinstance(a, (tuple, list)) and isinstance(b, (tuple, list)):
        arg_ = int((len(a) + len(b)) / 2)
        if arg_ == 3:
            res = [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                   a[0] * b[1] - a[1] * b[0]]
        elif arg_ == 2:
            res = a[0] * b[1] - b[0] * a[1]
        else:
            raise ValueError(
                "Vector must have a dimension of either 2 or 3."
            )

    return res


def deg_rad(
        angle,
):
    angle = valid_num(angle)

    return angle * math.pi / 180


def rad_deg(
        angle,
):
    angle = valid_num(angle)

    return angle * 180 / math.pi


def pow(
        x,
        n=2,
):
    x = valid_num(x)
    n = valid_num(n)

    if n == 0:
        return 1
    else:

        if x < 0:
            if n % 2 == 0:
                raise ValueError(
                    'Error: Only non-even roots can have a negative argument.'
                )
            else:
                return x ** n
        else:
            return x ** n


def quad_root(
        a,
        b,
        c,
):
    for arg_name in (a, b, c):
        if not isinstance(arg_name, (int, float)):
            raise TypeError(
                f"Argument '{arg_name}'must be an integer or float."
            )

    deter = pow(b ** 2 - 4 * a * c, 0.5)

    if a == 0:
        raise ValueError(
            'Error: Must be a quadratic to use dummy.'
        )
    else:

        x_0 = (-b + deter) / 2 * a
        x_1 = (-b - deter) / 2 * a

    return x_0, x_1


def is_even(
        number,
):
    number = valid_num(number)

    if number % 2 == 0:
        return True
    else:
        return False


def is_odd(
        number,
):
    if is_even(number) == False:
        return True
    else:
        return False


def dice_roll(
        dice_num=1,
):
    dice_num = valid_num(dice_num)
    res = []
    for i in range(dice_num):
        roll = random.randint(1, 6)
        res.append(roll)

    return res
