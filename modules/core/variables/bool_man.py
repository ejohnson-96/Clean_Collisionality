
def valid_bool(
        a,
):
    if not isinstance(a, bool):
        raise ValueError(
            f"Error: Argument {a} needs to be a boolean type, "
            f"instead type of {type(a)}."
        )
    else:
        return a


def and_(
        a,
        b,
):
    if valid_bool(a) == 1 and valid_bool(b) == 1:
        return True
    else:
        return False


def nand_(
        a,
        b,
):
    if valid_bool(a) == 1 and valid_bool(b) == 1:
        return False
    else:
        return True


def or_(
        a,
        b,
):
    if valid_bool(a) == 1 or valid_bool(b) == 1:
        return True
    else:
        return False


def xor_(
        a,
        b,
):
    if valid_bool(a) != valid_bool(b):
        return True
    else:
        return False


def not_(
        a,
):
    a = valid_bool(a)
    return not a


def nor_(
        a,
        b,
):
    a = valid_bool(a)
    b = valid_bool(b)

    if (a == 0) and (b == 0):
        return True
    elif (a == 0) and (b == 1):
        return False
    elif (a == 1) and (b == 0):
        return False
    elif (a == 1) and (b == 1):
        return True
    else:
        raise ValueError(
            f"Error: Logic cannot be determined, arguments passed "
            f"are {a} and {b} of types {type(a)} and {type(b)} "
            "respectively."
        )


def xnor_(
        a,
        b,
):
    if valid_bool(a) == valid_bool(b):
        return True
    else:
        return False
