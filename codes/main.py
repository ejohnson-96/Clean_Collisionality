from core.time import tictoc as stopwatch
from core.constants import initialise_constants
from constants import *

stopwatch.start_time()
initialise_constants()

p = 'proton'
a = 'alpha'

particle_list = [p,a]
valid_enc = [6,7]
print('Currently loaded encounters:', valid_enc, '\n')

valid_full = ['f','F','full','Full','ful','Ful', 'ff', 'FF', 'fF', 'Ff']
valid_single = ['s', 'S', 'single', 'Single', 'ss', 'SS', 'sS', 'Ss']

h = 1
while h >0:
    data_set_input = input('Full data set or singular? (F/S)')
    if data_set_input in valid_full:
        enc = 0
        h = 0
    elif data_set_input in valid_single:
        h = 0
        g = 1
        while g > 0:
            enc_input = input('Please enter an encounter:')

            if enc_input in valid_enc:
                enc = enc_input
                g = 0
            elif enc_input == '':
                print('Error: No input provided.')
            else:
                print('Error: Argument provided is not valid,'
                      f' argument {enc_input}, is of type {type(enc_input)}.')
    else:
        print('Error: Please make a valid selection.')

dum_val_ = file_dir_gen(enc, valid_enc)






print('test')




stopwatch.end_time()