def const(

):
    const.x_dim = 10
    const.y_dim = 10
    const.title_size = 30
    const.label_size = 30
    const.tick_size = 24
    const.legend_size = 30
    const.font_family = 'Courier New'
    const.line_width = 2

    smooth = 15
    const.smooth = int(smooth/3)
    const.bin_width = 0.2

    return


def system_const(

):

    system_const.alphabet = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
            'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x',
            'y', 'z']
    system_const.numerical = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '0']
    system_const.symbols = ['!', '@', '#', '$', '%', '^', '&', '*', '(', ')']
    system_const.nato = ['alpha', 'beta', 'charlie', 'delta', 'echo',
             'foxtrot', 'golf', 'hotel', 'india',
             'juliet', 'kilo', 'lima', 'mike', 'november',
             'oscar', 'papa', 'quebec', 'romeo', 'sierra',
             'tango', 'uniform', 'victor', 'whiskey', 'x-ray',
             'yankee', 'zulu']

    return


def initialise_constants(

):
    const()
    system_const()
    res = 'Note: Constant initialisation successful.'

    return res


