import numpy as np

def smooth(
        data,
        box_pts=1,
):

    if not isinstance(box_pts, int):
        raise TypeError(
            "Error: Argument 'box_pts' must be an integer, "
            f"instead got {box_pts} of type {type(box_pts)}."
        )

    box = np.ones(box_pts)/box_pts
    data_smooth = np.convolve(data, box, mode='same')

    return data_smooth