import numpy as np
import os
import matplotlib.pyplot as plt
# import seaborn as sns
from pandas import read_csv as read_csv

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
dirs = [os.path.join(source_dir,d) for d in ld if "runs" in d and os.path.isdir(os.path.join(source_dir,d))]
dirs.sort()

for c, cond in enumerate(dirs):

    #predictions = read_csv(os.path.join(cond,"predictions.txt"), sep=' ')
    parameters = read_csv(os.path.join(cond,"parameters.txt"), sep='|')
    measurements = read_csv(os.path.join(cond,"measurements.txt"), sep=' ')
    dist = manual_parsing(os.path.join(cond,"results.txt"), ' ',float)
    t = np.linspace(0,40,101)
    with open(os.path.join(cond,"x_var.txt")) as ins:
        var_ax, var_color = ins.readline().split()


    for f, f_par in enumerate(np.unique(parameters[var_color])):

        plt.figure(figsize=[6, 6])
        plt.title(var_color + " " + str(f_par))
        ind = parameters[var_color]==f_par
        set_color(parameters[var_ax][ind])
        plt.plot(t, (dist[ind,:].T-dist[ind,0]))
        plt.xlabel("Time (s)")
        plt.ylabel("Distance (um)")
        plt.xlim(xmin=0)
        plt.ylim(ymin=0)
        plt.savefig(os.path.join(cond, "dynamics_" + str(f) + ".svg"))

    plt.show()




