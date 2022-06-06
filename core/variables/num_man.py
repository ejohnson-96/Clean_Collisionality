import math


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
        if not isinstance(arg_name, (int, float)):
            raise TypeError(
                f"Argument '{arg_name}'must be an integer or float."
            )

    res = a + b + c + d + e + f + g + h + i + j + k + l + m

    return res


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
        if not isinstance(arg_name, (int, float)):
            raise TypeError(
                f"Argument '{arg_name}'must be an integer or float."
            )

    res = a - b - c - d - e - f - g - h - i - j - k - l - m

    return res


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
    for arg_name in (a, b):
        if not isinstance(arg_name, (int, float)):
            raise TypeError(
                f"Argument '{arg_name}'must be an integer or float."
            )

    if b == 0:
        raise ValueError(
            "Second argument cannot be zero dummy."
        )

    res = a / b

    return res


def dot_product(
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
    if not isinstance(angle, (int, float)):
        raise ValueError(
            "Angle must be either an integer or float, instead "
            f"got {type(angle)}."
        )

    arg_ = angle * math.pi / 180

    return arg_


def rad_deg(
        angle,
):
    if not isinstance(angle, (int, float)):
        raise ValueError(
            "Angle must be either an integer or float, instead "
            f"got {type(angle)}."
        )

    arg_ = angle * 180 / math.pi

    return arg_


def pow(
        x,
        n=2,
):
    for arg_name in (x, n):
        if not isinstance(arg_name, (int, float)):
            raise TypeError(
                f"Argument '{arg_name}'must be an integer or float."
            )

    if n == 0:
        res = 1
    else:

        if x < 0:
            if n % 2 == 0:
                raise ValueError(
                    'Error: Only non-even roots can have a negative argument.'
                )
            else:
                res = x ** n
        else:
            res = x ** n

        res = x ** n

    return res


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
