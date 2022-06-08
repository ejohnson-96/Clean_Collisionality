from modules.core.time import tictoc as stopwatch
from modules.core.constants import initialise_constants

from constants import *
from modules.core.system import input as inpt
from modules.collisions.loadsave import loadsave as rw
from modules.core.variables import num_man as nm

stopwatch.start_time()
initialise_constants()

t = 'time'
p = 'proton'
a = 'alpha'

particle_list = [p,a]
valid_enc = [4,6,7]
print('Currently loaded encounters:', valid_enc, '\n')

valid_full = ['f','F','full','Full','ful','Ful', 'ff', 'FF', 'fF', 'Ff']
valid_single = ['s', 'S', 'single', 'Single', 'ss', 'SS', 'sS', 'Ss']

h = 1
while h >0:
    data_set_input = input('Full data set or singular? (F/S)')
    if data_set_input in valid_full:
        encount = 0
        h = 0
    elif data_set_input in valid_single:
        h = 0
        g = 1
        while g > 0:
            enc_input = input('Please enter an encounter:')
            if inpt.validate_input_number(enc_input):
                if int(enc_input) in valid_enc:
                    encount = int(enc_input)
                    g = 0
                elif enc_input == '':
                    print('Error: No input provided.')
                else:
                    print('Error: No coresponding encounter available.')
            else:
                print('Error: Argument provided is not valid,'
                      f' argument {enc_input}, is of type {type(enc_input)}.')
    else:
        print('Error: Please make a valid selection.')

# Generate directory strings for encounters
enc(encount, valid_enc)
print('\n')
error_files = enc.error_files_loaded
mm_data = rw.encounter_import(encount, valid_enc, error_files)
if error_files:
    error_data = rw.error_import()
sc_data = rw.sc_import()
print('Data import sucessful.')

# Generate equal sized arays for all data

mm_len_max = 0
error_len_max = 0
sc_len_max = 0

for x in enc.encounter:
    for y in mm_data[x].keys():
        if len(mm_data[x][y][t]) > mm_len_max:
            mm_len_max = len(mm_data[x][y][t])
            arg_encount_ = x
            arg_mm_file = y

    if error_files:
        for y in error_data.keys():
            if len(error_data[x][y][t]) > error_len_max:
                error_len_max = len(error_data[x][y][t])
                arg_encount_ = x
                arg_error_file = y
    else:
        pass

    for y in sc_data[x].keys():
        if len(sc_data[x][y][t]) > sc_len_max:
            sc_len_max = len(sc_data[x][y][t]) #here
            arg_encount_ = x
            arg_sc_file = y

lengths = [mm_len_max, error_len_max, sc_len_max]
max_len = max(lengths)

if max_len == lengths[0]:
    t_ = mm_data[arg_encount_][arg_mm_file][t]
elif max_len == lengths[1]:
    t_ = error_data[arg_encount_][arg_error_file][t]
elif max_len == lengths[2]:
    t_ = sc_data[arg_encount_][arg_sc_file][t]

for x in enc.encounter:
    for y in mm_data[x].keys():
        xp = mm_data[x][y][t]
        for z in mm_data[x][y].keys():
            fp = mm_data[x][y][z]
            mm_data[x][y][z] = np.interp(t_, xp, fp)

    if error_files:
        for y in error_data[x].keys():
            xp = error_data[x][y][t]
            for z in error_data[x][y][z]:
                fp = error_data[x][y][z]
                error_data[x][y][z] = np.interp(t_, xp, fp)

    for y in sc_data[x].keys():
        xp = sc_data[x][y][t]
        for z in sc_data[x][y].keys():
            fp = sc_data[x][y][z]
            sc_data[x][y][z] = np.interp(t_, xp, fp)

# Generate file and/or combine files
print('Generating data file... \n')
solar_data = {}
errors = {}
spc_data = {}

for particle in particle_list:
    solar_data[particle] = {}
    errors[particle] = {}

for z in mm_data[encout][enc.encounter_names[0]].keys():
    solar_data[p][z] = []
for z in mm_data[encount][enc.encounter_names[1]].keys():
    solar_data[a][z] = []

if error_files:
    for z in error_data[encount][enc.encounter_names[0]].keys():
        errors[p][z] = []
    for z in error_data[encount][enc.encounter_names[1]].keys():
        errors[a][z] = []

for y in range(len(enc.sc_names)):
    spc_data[enc.sc_names[y]] = {}
    for z in sc_data[encount][enc.sc_names[y]].keys():
        spc_data[enc.enc.sc_names[y]][z] = []

for x in range(enc.num_of_encs):
    encout = enc.encounter[0]
    for particle in particle_list:
        for z in solar_data[particle].keys():
            for w in range(len(mm_data[encount][enc.encounter_names[2*x]][z])):
                if nm.is_even
                solar_data[particle][z].append(mm_data[encount][enc.encounter_names[2*x + arg_x_]][z][w])






stopwatch.end_time()



