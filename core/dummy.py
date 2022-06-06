import constants as const
from core.variables import char_man as cm
from core.variables import string_man as sm
from system import fonts as font
from core.time import tictoc as stopwatch

const.initialise_constants()
stopwatch.start_time()



x = font.list_fonts()
y = x['name'][10]
print(y)





test = 'fsfjflsj i09sdflkj'
y = cm.nato_alphabet(test, False)
print(y)

res = cm.remove_char(y, 'a')

print(res)

stopwatch.end_time()




