
from plot_settings import *
import os
from pandas import read_csv
import numpy as np

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

source_dir = "../Fig4A_parameters_diffusion/data_summary/runs"
ld = os.listdir(source_dir)


marker = new_fig()

parameters = read_csv(os.path.join(source_dir,"parameters.txt"), sep='|')
measurements = read_csv(os.path.join(source_dir,"measurements.txt"), sep='|')

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

plt.scatter(st_out, D_out,label="Best fit: simulations",marker=marker.next(),alpha=0.6)


plt.subplots_adjust(left=0.15,bottom=0.15)

plt.xlim(xmin=0)
plt.ylim([0,0.06])
plt.xlabel("Stiffness ($pN/\mu m$)")
plt.ylabel("$D_1$ ($\mu m^2/s$)")


plt.legend(frameon=False)
plt.savefig("fig4/fig4A.pdf")
plt.show()




