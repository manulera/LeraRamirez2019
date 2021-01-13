
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
gillespie = np.genfromtxt("gillespie_results.txt",dtype=float,delimiter=",")
plt.figure()
for c, cond in enumerate(dirs):

    #predictions = read_csv(os.path.join(cond,"predictions.txt"), sep=' ')
    parameters = read_csv(os.path.join(cond,"parameters.txt"), sep='|')

    dist = manual_parsing(os.path.join(cond,"results.txt"), ' ',float)
    dist = dist-parameters["length"][0]/2.



    t = np.linspace(0,100,101)
    with open(os.path.join(cond,"x_var.txt")) as ins:
        var_ax, var_color = ins.readline().split()

    stiffness_all = list()
    D_all = list()
    for f, f_par in enumerate(np.unique(parameters[var_color])):
        stiffness_all.append(f_par)


        ind = parameters[var_color]==f_par
        data = dist[ind,:].ravel()
        mu, std = norm.fit(data)


        # Simulation data
        x,y = ecdf(data)



        # Normal fit
        # x = np.linspace(-1,1)
        # y = norm.cdf(x, mu, std)
        # plt.plot(x,y)

        t = 5.
        D = np.mean(np.power(data,2))/2./t
        D_all.append(D)
        x2 = np.linspace(x[0],x[-1],1000)
        prob_x = 1./np.sqrt(4*np.pi*D*t)*np.exp(-np.power(x2,2)/(4*D*t))
        cumprob_x = np.cumsum(prob_x)/np.sum(prob_x)
        if False:
            plt.title(var_color + " " + str(f_par) + cond.split("/")[-1])
            plt.plot(x, y)
            plt.plot(x2,cumprob_x)


    plt.scatter(stiffness_all,D_all,label="Cytosim")

plt.xlabel("Stiffness of the linker")
plt.ylabel("Couple diffussion")
plt.plot(gillespie[0,:],gillespie[1,:],label="Gillespie solution")
plt.legend()


plt.ylim(ymin=0)
plt.show()




