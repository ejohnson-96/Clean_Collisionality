

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
