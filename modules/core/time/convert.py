import datetime


def epoch_time(
        epoch_time,
):
    if not isinstance(epoch_time, (int, float)):
        raise TypeError(
            "Error: Argument must be either a float or an"
            f" integer, instead got type {type(epoch_time)}."
        )
    res = datetime.datetime.fromtimestamp(epoch_time)

    return res


