import validators

def valid_string(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Error: Argument {string} must be of type string, "
            f"instead got type of {type(string)}"
        )
    else:
        return string

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
    a = valid_string(a)
    b = valid_string(b)
    z = ' '.join([a, b])

    for arg_name in (c, d, e, f, g, h, i, j, k, l, m):
        if not isinstance(arg_name, (str, bool)):
            raise TypeError(
                f"Error: Argument '{arg_name}' must be a string,"
                f" instead got type {type(arg_name)}."
            )

    for default_val in (c, d, e, f, g, h, i, j, k, l, m):
        if not default_val:
            break
        else:
            z = ' '.join([z, default_val])

    return z


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
    a = valid_string(a)
    b = valid_string(b)
    z = ''.join([a, b])

    for arg_name in (c, d, e, f, g, h, i, j, k, l, m):
        if not isinstance(arg_name, (str, bool)):
            raise TypeError(
                f"Error: Argument '{arg_name}' must be a string,"
                f" instead got type {type(arg_name)}."
            )

    for default_val in (c, d, e, f, g, h, i, j, k, l, m):
        if not default_val:
            break
        else:
            z = ''.join([z, default_val])

    return z


def slash_check(
        string,
):
    string = valid_string(string)

    if string.endswith("/") or string.endswith(str("\\")):
        return string
    else:
        if '/' in string:
            return jwos(string, '/')
        elif "\\" in string:
            return jwos(string, "\\")
        else:
            raise ValueError(
                "Error: No slash was found in this string: "
                f"{string}, please try again."
            )


def web_check(
        url,
):
    url = valid_string(url)

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

    if validators.url(url):
        return url
    else:
        raise ValueError(
            f"Error: An the URL provided, {url},"
            " is invalid, please try again."
        )


def split(
        string,
        character=" ",
):
    string = valid_string(string)
    character = valid_string(character)

    if len(character) > 1:
        raise ValueError(
            f"Error: The split character provided {character}, "
            f"must be of type string, instead got {type(character)}."
        )
    else:
        return string.split(character)
