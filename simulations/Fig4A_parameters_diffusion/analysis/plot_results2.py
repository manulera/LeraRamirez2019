
import numpy as np
import os
import matplotlib.pyplot as plt
# import seaborn as sns
from pandas import read_csv as read_csv
from scipy.stats import norm
# sns.set(font_scale=1.5)
# sns.color_palette("viridis")
# sns.set_style("ticks")


color = ['red', 'green', 'blue', 'black']
colormap = plt.cm.viridis

def sort_plot(x,y,**kwargs):
    x = np.array(x)
    y = np.array(y)
    ind = np.argsort(x)
    plt.plot(x[ind],y[ind],**kwargs)

def ecdf(d):
    x = np.sort(d)
    y = np.linspace(0,1,len(d))
    return x,y


def manual_parsing(filename, delim, dtype,first_line=False):
    out = list()
    lengths = list()
    column_names = None
    with open(filename, 'r') as ins:
        if first_line:
            line = ins.readline().strip()
            column_names =line.split(delim)
        for line in ins:
            l = line.split(delim)
            out.append(l)
            lengths.append(len(l))
    lim = np.max(lengths)
    for l in out:
        while len(l) < lim:
            l.append("nan")
    if not first_line:
        return np.array(out, dtype=dtype)
    else:
        return np.array(out, dtype=dtype), column_names

def set_color(data):
    mini = np.min(data)
    maxi = np.max(data)
    color = [colormap(i) for i in (data-mini) / (maxi-mini)]
    plt.gca().set_color_cycle(color)

source_dir = "../data_summary/"
ld = os.listdir(source_dir)
dirs = [os.path.join(source_dir,d) for d in ld if "run" in d and os.path.isdir(os.path.join(source_dir,d))]
dirs.sort()
# gillespie = np.genfromtxt("gillespie_results.txt",dtype=float,delimiter=",")

for c, cond in enumerate(dirs[-1:]):


    parameters = read_csv(os.path.join(cond,"parameters.txt"), sep='|')
    measurements = read_csv(os.path.join(cond,"measurements.txt"), sep='|')
    data = np.genfromtxt('lansky2015_Fig3B.csv', dtype=float, delimiter=",", skip_header=6)


    var_color = "stifness"
    plt.figure()
    for f, f_par in enumerate(np.unique(parameters[var_color])):


        ind = parameters[var_color]==f_par

        sort_plot(parameters["diffusion"][ind],measurements["couple_diff"][ind])








plt.show()




