from ..variables import string_man as sm
from ..variables import char_man as cm
from PIL import Image

convert_message = 'File conversion complete.'

def png_jpg(
        file_name,
        loc,
):

    if not isinstance(loc, str):
        raise TypeError(
            "Error: Path location must be a string."
        )

    if not isinstance(file_name, str):
        raise TypeError(
            "Error: File name must be a string."
        )

    parent_path = sm.slash_check(loc)
    if file_name.endswith(".png"):
        path = sm.jwos(parent_path, file_name)
    elif '.' in file_name:
        raise ValueError(
            "Error: Invalid file type."
        )
    else:
        arg_ = sm.jwos(file_name, ".png")
        path = sm.jwos(parent_path, arg_)


    im = Image.open(path)
    arg_ = path
    for i in range(4):
        arg_ = cm.remove_end(arg_)
    target_name = sm.jwos(arg_,".jpg")
    rgb_im = im.convert('RGB')
    rgb_im.save(target_name)

    return convert_message