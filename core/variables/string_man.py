
def jws(
        a,
        b,
        c=False,
        d=False,
        e=False,
        f=False,
        g=False,
        h=False,
        i=False,
        j=False,
        k=False,
        l=False,
        m=False,
):
    if not isinstance(a and b, str):
        raise TypeError(
            f"Arguments must be strings, instead got string_a: {a.type} "
            f"and string_b: {b.type}."
        )

    for arg_name in (c, d, e, f, g, h, i, j, k, l, m):
        if not isinstance(arg_name, (str, bool)):
            raise TypeError(
                f"Argument '{arg_name}' must be a string,"
                f" instead got type {arg_name.type}."
            )

    arg_ = ' '.join([a, b])

    for default_val in (c, d, e, f, g, h, i, j, k, l, m):
        if not default_val:
            break
        else:
            arg_ = ' '.join([arg_, default_val])

    res = arg_

    return res


def jwos(
        a,
        b,
        c=False,
        d=False,
        e=False,
        f=False,
        g=False,
        h=False,
        i=False,
        j=False,
        k=False,
        l=False,
        m=False,
):
    if not isinstance(a and b, str):
        raise TypeError(
            f"Arguments must be strings, instead got string_a: {a.type} "
            f"and string_b: {b.type}."
        )

    for arg_name in (c, d, e, f, g, h, i, j, k, l, m):
        if not isinstance(arg_name, (str, bool)):
            raise TypeError(
                f"Argument '{arg_name}' must be a string,"
                f" instead got type {arg_name.type}."
            )

    arg_ = ''.join([a, b])

    for default_val in (c, d, e, f, g, h, i, j, k, l, m):
        if not default_val:
            break
        else:
            arg_ = ''.join([arg_, default_val])

    res = arg_

    return res


def slash_check(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, "
            f"instead got {type(string)}"
        )


    if string.endswith("/") or string.endswith(str("\\")):
        res = string
    else:
        if '/' in string:
            res = jws(string, '/')
        elif "\\" in string:
            res = jwos(string,"\\")
        else:
            raise ValueError(
                "Error: No slash in string."
            )

    return res

def split(
        string,
        character,
        loc,
):

    return