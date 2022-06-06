def valid_bool(
        a,
        b,
):
    if isinstance(a, bool) and isinstance(b, bool):
        return a, b
    else:
        raise ValueError(
            "Error: Only boolean arguments are allowed."
        )


def _and(
        a,
        b,
):
    a, b = valid_bool(a, b)

    if a == 1 and b == 1:
        return True
    else:
        return False


def _nand(
        a,
        b,
):
    a, b = valid_bool(a, b)

    if a == 1 and b == 1:
        return False
    else:
        return True


def _or(
        a,
        b,
):
    a, b = valid_bool(a, b)

    if a == 1 or b == 1:
        return True
    else:
        return False

