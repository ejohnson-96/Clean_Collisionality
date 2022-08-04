import validators

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
            res = jwos(string, '/')
        elif "\\" in string:
            res = jwos(string,"\\")
        else:
            raise ValueError(
                "Error: No slash in string."
            )

    return res

def web_check(
        url,
):
    if not isinstance(url, str):
        raise TypeError(
            "Error: URL argument must be a string, instead "
            f"got type {type(url)}"
        )

    http = 'https://'
    www = 'www.'

    if url[:7] == http:
        pass
    elif url[:3] == 'w':
        url = http + url
    else:
        url = http + www + url

    if url[-1] != '/':
        url = url + '/'

    valid_url = validators.url(url)

    if valid_url:
        return url
    else:
        raise ValueError(
            "Error: An invalid URL was provided."
        )



def split(
        string,
        character,
        loc,
):

    return

