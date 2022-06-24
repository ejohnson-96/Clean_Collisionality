import numpy as np


def smooth(
        data,
        box_pts=1,
):
    if not isinstance(data, np.ndarray):
        raise TypeError(
            "Error: Argument 'data' must be a numpy array, "
            f"instead got type of {type(data)}."
        )

    if not isinstance(box_pts, int):
        raise TypeError(
            "Error: Argument 'box_pts' must be an integer, "
            f"instead got {box_pts} of type {type(box_pts)}."
        )

    if box_pts == 0:
        return data
    else:
        box = np.ones(box_pts) / abs(box_pts)
        return np.convolve(data, box, mode='same')
