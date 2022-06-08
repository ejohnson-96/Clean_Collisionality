from modules.core.time import tictoc as stopwatch
from modules.core.constants import initialise_constants

from constants import *
from modules.core.system import input as inpt
from modules.collisions.loadsave import loadsave as rw

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

mm_data = rw.encounter_import(encount, valid_enc)
error_data = rw.error_import()
sc_data = rw.sc_import()






print('test')




stopwatch.end_time()



