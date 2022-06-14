

def validate_input_number(
        input,
):
    if not isinstance(input, str):
        raise TypeError(
            "Error: User input must a string, instead "
            f"got {input} of type {type(input)}."
        )
    try:
        # Convert it into integer
        val = int(input)
        return True
    except ValueError:
        try:
            # Convert it into float
            val = float(input)
            return True
        except ValueError:
            return False


def return_all_keys(
        dictionary,
):
    if not isinstance(dictionary, dict):
        raise TypeError(
            'Error: Argument must be a dictionary, '
            f'instead got type {type(dictionary)}'
        )
    result = []
    for key, value in dictionary.items():
        if type(value) is dict:
            new_keys = return_all_keys(value)
            result.append(key)
            for innerkey in new_keys:
                result.append(f'{key}/{innerkey}')
        else:
            result.append(key)
    return result
