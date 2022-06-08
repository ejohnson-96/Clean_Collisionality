from modules.core.variables import string_man as sm
from modules.core.constants import system_const

system_const()


def capital_first_letter(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, got {type(string)}"
        )
    result = string.capitalize()

    return result


def capital_all_first_letter(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, got {type(string)}"
        )

    result = ' '.join(elem.capitalize() for elem in string.split())

    return result


def capital_all_letter(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, got {type(string)}"
        )

    result = string.upper()

    return result


def lower_first_letter(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, got {type(string)}"
        )

    result = string[0].lower() + string[1:]

    return result


def lower_all_first_letter(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, got {type(string)}"
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
            f" got {type(string)}."
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
            f" got {type(string)}."
        )

    result = string[1:]

    return result


def replace_char(
        string,
        character,
        replacement,
):
    arg_name_ = [string, character, replacement]
    for arg_ in arg_name_:
        if not isinstance(arg_, str):
            raise ValueError(
                f"Argument {arg_}, is not of string type, "
                f"instead got type {type(arg_)}."
            )

    if character in string:
        res = string.replace(character, replacement)
    else:
        raise TypeError(
            f"Character {character}, does not appear in the "
            f"string provided, {string}."
        )

    return res


def remove_char(
        string,
        character,
):
    res = replace_char(string, character, '')

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


def nato_alphabet(
        string,
        _list=False,
):
    alphabet = system_const.alphabet
    numerical = system_const.numerical
    symbols = system_const.symbols
    nato_alphabet = system_const.nato

    if not isinstance(string, str):
        raise ValueError(
            "Argument passed is a not string, instead "
            f" got {string.type}."
        )

    phon_test = string.split()

    if lower_all_letter(phon_test[0]) in nato_alphabet:
        phonetic = True
    else:
        phonetic = False

    res = ''
    res_list = []

    prev = None

    if not phonetic:
        first = True
        for letter in string:
            if lower_all_letter(letter) in alphabet:
                arg_ = alphabet.index(lower_all_letter(letter))
                if first:
                    res = sm.jwos(res, nato_alphabet[arg_])
                    res_list.append(capital_first_letter(nato_alphabet[arg_]))
                else:
                    res = sm.jws(res, nato_alphabet[arg_])
                    res_list.append(capital_first_letter(nato_alphabet[arg_]))
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
            if lower_all_letter(char) in nato_alphabet:
                arg_ = nato_alphabet.index(lower_all_letter(char))
                if first:
                    res = sm.jwos(res, alphabet[arg_])
                    res_list.append(
                        capital_first_letter(capital_first_letter(alphabet[arg_])))
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


def validate_case_entry(
        string,
        single=False,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Argument passed is not a string, got {type(string)}"
        )

    if single:
        if len(string) > 1:
            raise ValueError(
                "Error: Only single characters can be passed as an "
                f"argument, instead got {string} of type {type(string)}."
            )
    else:
        pass

    return string


def is_upper(
        char,
        single=False,
):
    validate_case_entry(char, single)
    res = char.isupper()

    return res


def is_lower(
        char,
        single=False,
):
    validate_case_entry(char, single)
    res = char.islower()

    return res

