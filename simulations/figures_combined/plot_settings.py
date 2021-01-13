
import matplotlib.pyplot as plt
from matplotlib import rcParams, cycler

import itertools

marker = None
# Color-blind friendly color cycle
CB_color_cycle = ['006BA4', 'FF800E', 'ABABAB', '595959','5F9ED1', 'C85200', '898989', 'A2C8EC', 'FFBC79', 'CFCFCF']


what = "article"

# rcParams["legend.frameon"]=False
rcParams['axes.prop_cycle'] = cycler(color=CB_color_cycle)

if what=="article":

    plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
    rcParams['lines.linewidth'] = 2
    rcParams['legend.fontsize']=12

    rcParams['legend.edgecolor']="black"
    rcParams["lines.markersize"]=8



elif what=="poster":

    plt.rc('axes', labelsize=22)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
    rcParams['lines.linewidth'] = 4
    rcParams['legend.fontsize']=16
    rcParams["lines.markersize"]=10

def new_fig():
    marker = itertools.cycle(('o', 's', '^'))

    if what=="article":
        plt.figure(figsize=[6, 6])
    elif what == "poster":
        plt.figure(figsize=[7, 7])

    return marker


def table_legend(handles,row_names, column_names,title_space=5,row_space=3,title_prespace=20,title_after_space=5,**kwargs):

    total_list =list()
    for l in handles:
        for ll in l:
            if type(ll)==list:
                total_list.append(ll[0])
            else:
                total_list.append(ll)

    title_space = " "*title_space

    title = title_prespace*" "+title_space.join(column_names) + title_after_space*" "

    row_names = [c + row_space*" " for c in row_names]

    leg = plt.legend(total_list, row_names + [''] * (len(total_list) - 2),
               title=title,
               ncol=len(column_names)-1, numpoints=1, markerfirst=False,**kwargs)

    plt.setp(leg.get_title(), fontsize=12)