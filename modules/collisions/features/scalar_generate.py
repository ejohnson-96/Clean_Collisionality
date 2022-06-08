import numpy as np

def scalar_velocity(
        data,
):
    t = 'time'
    p = 'proton'
    a = 'alpha'

    data[p]['v_mag'] = []
    data[a]['v_mag'] = []

    L_p = len(data[p][t])
    L_a = len(data[a][t])

    for i in range(L_p):
        val = (data[p]['vp1_x'][i]) ** 2 + (data[p]['vp1_y'][i]) ** 2 + (data[p]['vp1_z'][i]) ** 2
        if val < 0:
            val = 0
        else:
            pass
        data[p]['v_mag'].append(np.sqrt(val))
    for i in range(L_a):
        val = (data[a]['va_x'][i]) ** 2 + (data[a]['va_y'][i]) ** 2 + (data[a]['va_z'][i]) ** 2
        if val < 0:
            val = 0
        else:
            pass
        data[a]['v_mag'].append(np.sqrt(val))

    return data