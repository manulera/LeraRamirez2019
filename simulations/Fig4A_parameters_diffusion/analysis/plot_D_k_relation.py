
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


gillespie = np.genfromtxt("results_D_stiffness_results.txt",float,delimiter=",")

for c, cond in enumerate(dirs[-1:]):


    parameters = read_csv(os.path.join(cond,"parameters.txt"), sep='|')
    measurements = read_csv(os.path.join(cond,"measurements.txt"), sep='|')

    stiffness = parameters["stifness"]
    D = measurements["couple_diff"]
    D0 = parameters["diffusion"]


    # Plot cytosim

    st_out = list()
    D_out = list()

    for st in np.unique(stiffness):
        st_out.append(st)
        ind = np.equal(stiffness, st)
        this_D = D[ind]
        i = np.argmin(np.abs(this_D - 0.011))

        D_out.append(D0[ind][i])

    plt.figure()

    plt.scatter(st_out, D_out,label="Best fit: Cytosim")

    st_out = list()
    D_out = list()
    stiffness, D, D0 = gillespie

    for st in np.unique(stiffness):
        st_out.append(st)
        ind = np.equal(stiffness, st)
        this_D = D[ind]
        i = np.argmin(np.abs(this_D - 0.011))

        D_out.append(D0[ind][i])

    plt.plot(st_out, D_out,label="Best fit: Gillespie")

plt.xlabel("Stiffness ($pN/\mu m$)")
plt.ylabel("Diffusion ($\mu m^2/s$)")


plt.legend()

plt.show()




