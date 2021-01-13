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

def sort_plot(x,y,fun,**kwargs):
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
    print cond
    predictions = read_csv(os.path.join(cond,"predictions.txt"), sep=' ')
    parameters = read_csv(os.path.join(cond,"parameters.txt"), sep='|')
    measurements = read_csv(os.path.join(cond,"measurements.txt"), sep=' ')

    with open(os.path.join(cond,"x_var.txt")) as ins:
        var_ax, var_color = ins.readline().split()


    plt.figure(figsize=[6, 6])
    for f, f_par in enumerate(np.unique(parameters[var_color])):

        ind = parameters[var_color] == f_par
        if cond!="../data_summary/runs4":
            sc = plt.scatter(parameters[var_ax][ind], measurements["speed"][ind],label=var_color + " " +str(f_par))
            sort_plot(parameters[var_ax][ind], predictions["speed"][ind],plt.plot,color=sc.get_facecolors()[0])
        else:
            sc = plt.scatter(parameters[var_ax][ind], measurements["speed"][ind], label=var_color + " " + str(f_par))
            sort_plot(parameters[var_ax][ind], predictions["speed"][ind], plt.semilogx, color=sc.get_facecolors()[0])

    plt.xlabel(var_ax)
    plt.ylabel("speed")
    plt.legend()
    plt.xlim(xmin=0)
    # plt.ylim([0.,0.05])
    plt.ylim(ymin=0)
    plt.savefig(os.path.join(cond,"plot.svg"))
    plt.title(cond)
plt.show()




