import numpy as np
import os
import matplotlib.pyplot as plt
# import seaborn as sns
from pandas import read_csv as read_csv

# sns.set(font_scale=1.5)
# sns.color_palette("viridis")
# sns.set_style("ticks")


color = ['red', 'green', 'blue', 'black']
colormap = plt.cm.jet

def sort_plot(x,y,fun=plt.plot,**kwargs):
    x = np.array(x)
    y = np.array(y)
    ind = np.argsort(x)
    fun(x[ind],y[ind],**kwargs)

source_dir = "../data_summary/"
ld = os.listdir(source_dir)
dirs = [os.path.join(source_dir,d) for d in ld if "run" in d and os.path.isdir(os.path.join(source_dir,d))]
dirs.sort()
#dirs = [os.path.join(source_dir,"runs1"),os.path.join(source_dir,"runs2")]
for c, cond in enumerate(dirs):

    predictions = read_csv(os.path.join(cond,"predictions.txt"), sep=' ')
    parameters = read_csv(os.path.join(cond,"parameters.txt"), sep='|')
    measurements = read_csv(os.path.join(cond,"measurements.txt"), sep=' ')

    with open(os.path.join(cond,"x_var.txt")) as ins:
        var_ax, var_color = ins.readline().split()

    plt.figure(figsize=[6,6])
    for f, f_par in enumerate(np.unique(parameters[var_color])):
        ind = parameters[var_color]==f_par

        sc = plt.scatter(parameters[var_ax][ind], measurements["bound_mts"][ind], label=var_color + " " + str(f_par))
        sort_plot(parameters[var_ax][ind], predictions["bound_mts"][ind],plt.semilogx, color=sc.get_facecolors()[0])
    plt.xlabel(var_ax)
    plt.ylabel("bound_mts")


    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.xlabel(var_ax)
    plt.ylabel("bound_mts")
    plt.legend()
    plt.savefig(os.path.join(cond, "plot_bound.svg"))
    plt.figure(figsize=[6, 6])

    for f, f_par in enumerate(np.unique(parameters[var_color])):

        ind = parameters[var_color] == f_par
        # sc = plt.scatter(parameters[var_ax][ind], measurements["force"][ind],label=var_color + " " +str(f_par))
        # sort_plot(parameters[var_ax][ind], predictions["force"][ind],color=sc.get_facecolors()[0])

        sc = plt.scatter(parameters[var_ax][ind], measurements["force"][ind], label=var_color + " " + str(f_par))
        sort_plot(parameters[var_ax][ind], 2*predictions["force"][ind],plt.semilogx, color=sc.get_facecolors()[0])
        sort_plot(parameters[var_ax][ind], 2 * predictions["force2"][ind], plt.semilogx, color=sc.get_facecolors()[0])
        # sort_plot(parameters[var_ax][ind], predictions["max_force"][ind], plt.semilogx, color=sc.get_facecolors()[0],linestyle="--")

    plt.xlabel(var_ax)
    plt.ylabel("force")
    plt.legend()
    plt.xlim(xmin=0)
    # plt.ylim([0.,0.05])
    plt.ylim(ymin=0)
    plt.title(cond)
    plt.savefig(os.path.join(cond,"plot.svg"))

    plt.show()




