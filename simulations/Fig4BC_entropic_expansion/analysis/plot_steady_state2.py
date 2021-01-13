import numpy as np
import os
import matplotlib.pyplot as plt
# import seaborn as sns
from pandas import read_csv as read_csv

# sns.set(font_scale=1.5)
# sns.color_palette("viridis")
# sns.set_style("ticks")
kT = 0.0042
def drag(D,c,L,a,visc):
    rho_c = c*a/L

    drag_MT = 3 * L / (np.log(0.5 * L / 0.0125) + 0.312) * np.pi * visc
    return c*kT/D/(1.-rho_c)/2.+drag_MT

def entropic(c,L,a):
    l = L/a

    return -kT / a * np.log((l + 1 - c) / (l + 1))

def speed_entropic(L,c,D,visc,a):
    ent = entropic(c,L,a)
    d = drag(D,c,L,a,visc)
    return ent/d



color = ['red', 'green', 'blue', 'black']
colormap = plt.cm.jet

def sort_plot(x,y,fun,**kwargs):
    x = np.array(x)
    y = np.array(y)
    ind = np.argsort(x)
    fun(x[ind],y[ind],**kwargs)

source_dir = "../data_summary/"
ld = os.listdir(source_dir)
dirs = [os.path.join(source_dir,d) for d in ld if "predicted_pairs" in d and os.path.isdir(os.path.join(source_dir,d))]
dirs.sort()
data = np.genfromtxt('lansky2015_Fig3B.csv',dtype=float,delimiter=",",skip_header=6)

dirs = [os.path.join(source_dir,'runs_chosen_nu4')]
#dirs = [os.path.join(source_dir,"runs1"),os.path.join(source_dir,"runs2")]
for c, cond in enumerate(dirs):
    print cond
    # predictions = read_csv(os.path.join(cond,"predictions.txt"), sep=' ')
    parameters = read_csv(os.path.join(cond,"parameters.txt"), sep='|')
    measurements = read_csv(os.path.join(cond,"measurements.txt"), sep=' ')

    with open(os.path.join(cond,"x_var.txt")) as ins:
        var_ax, var_color = ins.readline().split()


    for f, f_par in enumerate(np.unique(parameters[var_color])):
        plt.figure(figsize=[6, 6])
        plt.scatter(data[:, 0], data[:, 1] / 1000, label="Lansky points")
        plt.title(f_par)
        ind = parameters[var_color] == f_par

        occup = np.array(parameters["occup"])
        ind = np.logical_and(ind,np.equal(occup,0.7))
        # print np.sum(ind)
        visc = np.array(parameters["viscosity"][ind])[0]
        sc = plt.scatter(parameters[var_ax][ind], measurements["speed"][ind],label=var_color + " " +str(f_par))
        # sc = plt.scatter(parameters[var_ax][ind], measurements["speed2"][ind], label=var_color + " " + str(f_par))
        # sc = plt.scatter(parameters[var_ax][ind], measurements["speed3"][ind], label=var_color + " " + str(f_par))
        L = np.linspace(0,10)

        # sp = speed_entropic(L, L*occup/0.008, 0.085,visc,0.008)
        # plt.plot(L, sp,c='black')

        plt.ylim(ymax=0.04)
        plt.xlim([0,10])
        # sort_plot(parameters[var_ax][ind], predictions["speed"][ind],plt.plot,color=sc.get_facecolors()[0])


        #sort_plot(parameters[var_ax][ind], predictions["speed"][ind], plt.semilogx, color=sc.get_facecolors()[0])

    plt.xlabel(var_ax)
    plt.ylabel("speed")
    plt.legend()
    # plt.xlim(xmin=0)
    # plt.ylim([0.,0.05])
    # plt.ylim(ymin=0)
    plt.savefig(os.path.join(cond,"plot.svg"))

plt.show()




