import math
from modules.core.features import smooth as smoothing

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from modules.core.loadsave import file_dir as fd
from modules.core.features import graph as graph
from modules.core.variables import char_man as cm, string_man as sm
import generate_files as gen_files

#vals = [0.1]
#for val in vals:
#gen_files.encounter_generator(val, False)

slash = fd.slash()
parent_name = fd.dir_name(fd.dir_parent())
print(parent_name)
path = fd.dir_parent()
for i in range(len(parent_name) + 1):
    path = cm.remove_end(path)
path = sm.jwos(sm.slash_check(sm.slash_check(path)),'data',slash,'save',slash)

encounters = ['EA', 'E4', 'E6', 'E7']
vals = [0.1,0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0]
paths = {}
for encount in encounters:
    paths[encount] = {}
    for val in vals:
        paths[encount][val] = sm.jwos(path,encount,slash,str(val),slash)

h = 1
user_enc = 'EA'
user_rad = 1.0

path1 = paths[user_enc][0.1]
path2 = paths['E6'][0.1]
path3 = paths[user_enc][1.0]

theta_i = np.loadtxt(sm.jwos(path3,'theta_f.txt'))
theta_f = np.loadtxt(sm.jwos(path2,'wind_theta.txt'))
theta_w = np.loadtxt(sm.jwos(path3,'theta_f.txt'))

data = [theta_i, theta_f, theta_w]

for file in data:
    for i in range(len(file)):
        if file[i] > 15:
            file[i] = 15
        else:
            file[i] = file[i]


label_i = '0.1 - 0.2 AU'
label_f = str(val) + ' AU'

X = np.linspace(0, 2, 1000)

def maxwellian(x, r, m,s):
    return (r / (s * np.sqrt(math.pi))) * np.exp(- (x - m) ** 2 / (2 * (s ** 2)))


def fit_function(x, r, m, s, q, n, t, ):
    return maxwellian(x, r, m,s) - maxwellian(x,q,n,t)

def second_fit_function(x, r, m, s, q, n, t, ):
    return maxwellian(x, r, m,s) + maxwellian(x,q,n,t)

val = 0.48
val2 = val
val3 = val

n = 1

bn_i = 100 #int((max(theta_i) - min(theta_i))/val)
bn_f = bn_i #int((max(theta_f) - min(theta_f))/val2)
bn_w = bn_i #int((max(theta_w)-min(theta_w))/val3)

theta_i = theta_i[np.logical_not(np.isnan(theta_i))]
theta_f = theta_f[np.logical_not(np.isnan(theta_f))]
theta_w = theta_w[np.logical_not(np.isnan(theta_w))]


print(bn_i, bn_f, bn_w)


for i in range(len(theta_i)):
    if not isinstance(theta_i[i], (float, int)):
        print(theta_i[i])

smooth_val = 1

weg1 = np.ones_like(theta_i)/float(len(theta_i))
weg2 = np.ones_like(theta_f)/float(len(theta_f))
weg3 = np.ones_like(theta_w)/float(len(theta_w))

data_entries_i, bins_i= np.histogram(smoothing.smooth(theta_i, smooth_val), range=(0,15), weights=weg1, bins=bn_i, density=1)
data_entries_f, bins_f = np.histogram(smoothing.smooth(theta_f,smooth_val), range=(0,15), weights=weg2, bins=bn_f, density=1)
data_entries_w, bins_w = np.histogram(smoothing.smooth(theta_w, smooth_val), range=(0,15), weights=weg3, bins=bn_w, density=1)

binscenters_i = np.array([0.5 * (bins_i[i] + bins_i[i + 1]) for i in range(len(bins_i) - 1)])
binscenters_f = np.array([0.5 * (bins_f[i] + bins_f[i + 1]) for i in range(len(bins_f) - 1)])
binscenters_w = np.array([0.5 * (bins_w[i] + bins_w[i + 1]) for i in range(len(bins_w) - 1)])





#pop_i, pcov_i = curve_fit(maxwellian, xdata=binscenters_i, ydata=data_entries_i, p0=[0.1,6,0.1])
#popt_f, pcov_f = curve_fit(fit_function, xdata=binscenters_f, ydata=data_entries_f,maxfev=80000, p0=[0.1,1.2,0.01,1,6,1])
#pop_w, pcov_w = curve_fit(second_fit_function, xdata=binscenters_w, ydata=data_entries_w, maxfev=180000,p0=[0.1,1.0,0.1,0.1,6,0.1] )

#xspace = np.linspace(0, 15, 1000)
#yspace = fit_function(xspace, *popt_f)
#zspace = maxwellian(xspace, *pop_i)
#wspace = second_fit_function(xspace, *pop_w)

x = np.linspace(0,14, bn_i)
y = {'1.0 AU':data_entries_i, 'Wind @ 1.0 AU':data_entries_f, }#'1.0 AU':data_entries_w}


y_label = 'Probability Density'
x_label = r'$\alpha$-Proton Relative Temperature'

color = ['gray','black',]#'blue']
style = ['-','--']#,'--']

graph.graph(x, y, colours=color, style_line=style, title='', x_axis=x_label, y_axis=y_label, limits=True, x_lim=14, y_lim=1, line_width=2)
