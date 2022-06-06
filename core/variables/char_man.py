from core.variables import string_man as sm


def capital_first_letter(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, got {string.type}"
        )
    result = string.capitalize()

    return result


def capital_all_first_letter(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, got {string.type}"
        )

    result = ' '.join(elem.capitalize() for elem in string.split())

    return result


def capital_all_letter(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, got {string.type}"
        )

    result = string.upper()

    return result


def lower_first_letter(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, got {string.type}"
        )

    result = string[0].lower() + string[1:]

    return result


def lower_all_first_letter(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, got {string.type}"
        )

    result = ' '.join(elem.lower() for elem in string.split())

    return result


def lower_all_letter(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, got {string.type}"
        )

    result = string.lower()

    return result


def remove_end(
        string,
):
    if not isinstance(string, str):
        raise ValueError(
            "Argument passed is a not string, instead "
            f" got {string.type}."
        )
    l = len(string)
    result = string[:l - 1]

    return result


def remove_begin(
        string,
):
    if not isinstance(string, str):
        raise ValueError(
            "Argument passed is a not string, instead "
            f" got {string.type}."
        )

    result = string[1:]

    return result

def remove_char(
        string,
        character,
):
    arg_name_ = []
    for arg_ in arg_name_:
        if not isinstance(arg_, str):
            raise ValueError(
                f"Argument {arg_}, is not of string type, "
                f"instead got type {type(arg_)}."
            )

    if character in string:
        res = string.replace(character, '')
    else:
        raise TypeError(
            f"Character {character}, does not appear in the "
            f"string provided, {string}."
        )

    return res


def add_space_front(
        string,
):
    if not isinstance(string, str):
        raise ValueError(
            "Argument passed is a not string, instead "
            f" got {string.type}."
        )

    res = string + ' '

    return res


def add_space_end(
        string,
):
    if not isinstance(string, str):
        raise ValueError(
            "Argument passed is a not string, instead "
            f" got {string.type}."
        )

    res = ' ' + string

    return res


alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
            'm', 'n', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x',
            'y', 'z']
numerical = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
symbols = ['!', '@', '#', '$', '%', '^', '&', '*', '(', ')']
nato_alph = ['alpha', 'beta', 'charlie', 'delta', 'echo',
             'foxtrot', 'golf', 'hotel', 'india',
             'juliet', 'kilo', 'lima', 'mike', 'november',
             'oscar', 'papa', 'quebec', 'romeo', 'sierra',
             'tango', 'uniform', 'victor', 'whiskey', 'x-ray',
             'yankee', 'zulu']


def nato_alphabet(
        string,
        _list=False,
):
    if not isinstance(string, str):
        raise ValueError(
            "Argument passed is a not string, instead "
            f" got {string.type}."
        )

    phon_test = string.split()
    if phon_test[0] in nato_alph:
        phonetic = True
    else:
        phonetic = False

    res = ''
    res_list = []

    prev = None

    if not phonetic:
        first = True
        for letter in string:
            if letter in alphabet:
                arg_ = alphabet.index(letter)
                if first:
                    res = sm.jwos(res, nato_alph[arg_])
                    res_list.append(capital_first_letter(nato_alph[arg_]))
                else:
                    res = sm.jws(res, nato_alph[arg_])
                    res_list.append(capital_first_letter(nato_alph[arg_]))
            elif letter in numerical:
                if first:
                    res = sm.jwos(res, letter)
                    res_list.append(letter)
                else:
                    if prev in numerical:
                        res = sm.jwos(res, letter)
                        res_list.append(letter)
                    else:
                        res = sm.jws(res, letter)
                        res_list.append(letter)
            elif letter in symbols:
                res = sm.jwos(res, letter)
                res_list.append(letter)
            elif letter == ' ':
                if prev in symbols:
                    res = res
                    res_list = res_list
                else:
                    res = res + ','
                    res_list.append(' ')
            else:
                raise ValueError(
                    'Error: Letter not in the alphabet,'
                    f' instead got {letter}.'
                )
            prev = letter
            first = False
    else:

        chars = string.split(' ')

        first = True
        prev = None
        for char in chars:
            if char in nato_alph:
                arg_ = nato_alph.index(char)
                if first:
                    res = sm.jwos(res, alphabet[arg_])
                    res_list.append(capital_first_letter(capital_first_letter(alphabet[arg_])))
                else:
                    res = sm.jws(res, alphabet[arg_])
                    res_list.append(' ')
                    res_list.append(capital_first_letter(alphabet[arg_]))
            else:
                second = True
                for s in char:
                    if s in numerical:
                        if first:
                            res = sm.jwos(res, s)
                            res_list.append(capital_first_letter(s))
                        else:
                            if prev in numerical:
                                if second:
                                    res = sm.jws(res, s)
                                    res_list.append(' ')
                                    res_list.append(s)
                                else:
                                    res = sm.jwos(res, s)
                                    res_list.append(capital_first_letter(s))
                            else:
                                res = sm.jws(res, s)
                                res_list.append(' ')
                                res_list.append(capital_first_letter(s))
                    elif s in symbols:
                        res = sm.jwos(res, s)
                        res_list.append(s)
                    else:
                        raise ValueError(
                            'Error: Letter not in the NATO alphabet,'
                            f' instead got {s}.'
                        )
                    second = False
                    prev = s
            first = False

    if _list:
        res_arg_ = res_list
    else:
        res_arg_ = capital_all_first_letter(res)

    return res_arg_

