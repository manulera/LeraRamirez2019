import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import rcParams, cycler
from pandas import read_csv as read_csv

import itertools

marker = None
# Color-blind friendly color cycle
CB_color_cycle = ['006BA4', 'FF800E', 'ABABAB', '595959','5F9ED1', 'C85200', '898989', 'A2C8EC', 'FFBC79', 'CFCFCF']



# plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('axes', labelsize=22)    # fontsize of the x and y labels
# plt.rc('xtick', labelsize=14)    # fontsize of the tick labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
# plt.rc('ytick', labelsize=14)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
rcParams['axes.prop_cycle']=cycler(color=CB_color_cycle)
# rcParams['lines.linewidth'] = 2
# rcParams['lines.linewidth'] = 2
rcParams['lines.linewidth'] = 4
# rcParams['legend.fontsize']=12
rcParams['legend.fontsize']=18
rcParams["legend.frameon"]=False
# rcParams["lines.markersize"]=8
rcParams["lines.markersize"]=10

def new_fig():
    marker = itertools.cycle(('o', 's', '^'))
    # plt.figure(figsize=[6, 6])
    plt.figure(figsize=[7, 7])

    return marker