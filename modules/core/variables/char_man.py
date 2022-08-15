from modules.core.variables import string_man as sm
from modules.core.constants import system_const

system_const()


def valid_char(
        string,
):
    if not isinstance(string, str):
        raise TypeError(
            f"Error: Argument {string} must be of type string, "
            f"instead got type of {type(string)}"
        )
    else:
        return string


def capital_first_letter(
        string,
):
    string = valid_char(string)
    return string.capitalize()


def capital_all_first_letter(
        string,
):
    string = valid_char(string)
    result = ' '.join(elem.capitalize() for elem in string.split())

    return result


def capital_all_letter(
        string,
):
    string = valid_char(string)

    return string.upper()


def lower_first_letter(
        string,
):
    string = valid_char(string)

    return string[0].lower() + string[1:]


def lower_all_first_letter(
        string,
):
    string = valid_char(string)
    result = ' '.join(elem.lower() for elem in string.split())

    return result


def lower_all_letter(
        string,
):
    string = valid_char(string)

    return string.lower()


def remove_end(
        string,
):
    string = valid_char(string)

    return string[:len(string) - 1]


def remove_begin(
        string,
):
    string = valid_char(string)

    return string[1:]


def replace_char(
        string,
        character,
        replacement,
):
    string = valid_char(string)
    character = valid_char(character)
    replacement = valid_char(replacement)

    if character in string:
        return string.replace(character, replacement)
    else:
        raise ValueError(
            f"Error: The character provided, {character}, does not "
            f"appear in the provided string, {string}."
        )


def remove_char(
        string,
        character,
):
    return replace_char(string, character, '')


def add_space_front(
        string,
):
    string = valid_char(string)

    return string + ' '


def add_space_end(
        string,
):
    string = valid_char(string)

    return ' ' + string


def nato_alphabet(
        string,
        _list=False,
):
    alphabet = system_const.alphabet
    numerical = system_const.numerical
    symbols = system_const.symbols
    nato_alphabet = system_const.nato

    string = valid_char(string)

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
        return res_list
    else:
        return capital_all_first_letter(res)


def validate_case_entry(
        string,
        single=False,
):
    string = valid_char(string)

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
