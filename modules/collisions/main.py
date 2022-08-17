from graphs import fig1, fig2, fig3, fig4, fig5
from modules.core.loadsave import file_dir as fd
from modules.core.variables import char_man as cm, string_man as sm
import generate_files as gen_files

# vals = [1.0]
# for val in vals:
#    gen_files.encounter_generator(val, False)

slash = fd.slash()


def gen_load_paths(

):
    parent_name = fd.dir_name(fd.dir_parent())
    path = fd.dir_parent()
    print(path)
    for i in range(len(parent_name) + 1):
        path = cm.remove_end(path)
    path = sm.jwos(sm.slash_check(path), 'data', slash, 'save', slash)

    encounters = ['EA', 'E4', 'E6', 'E7']
    vals = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0]
    paths = {}
    for encount in encounters:
        paths[encount] = {}
        for val in vals:
            paths[encount][val] = sm.jwos(path, encount, slash, str(val), slash)

    return paths


paths = gen_load_paths()

enc = 'EA'
bin_num = 150
test = (bin_num/2)*10**-3
value = 0.03

fig1.fig1(enc, bin_num, value)
fig2.fig2(enc, bin_num, value)
fig3.fig3(enc, bin_num, value)
fig4.fig4(enc, bin_num, value)
fig5.fig5()
