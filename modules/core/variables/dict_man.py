
def valid_dict(
        input,
):
    if not isinstance(input, dict):
        raise TypeError(
            "Error: Arugment provided must be of type dict, "
            f"instead got type of {type(input)}."
        )
    else:
        return input


def recursive_items(
        dictionary
):
    for key, value in dictionary.items():
        if type(value) is dict:
            yield from recursive_items(value)
        else:
            yield key, value


def all_keys(
        dictionary,
):
    dictionary = valid_dict(dictionary)
    keys = []
    values = []
    for key, value in recursive_items(dictionary):
        keys.append(key)
        values.append(value)

    return keys, values


def sort_alpha(
        dictionary,
):
    dictionary = valid_dict(dictionary)
    return {key: value for key, value in sorted(dictionary.items())}


