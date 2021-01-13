
#!/usr/bin/env python

import sys, os
import numpy as np
from pandas import read_csv

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

def calc_speed(t,y):

    dims = y.shape
    out = list()
    for i in range(dims[0]):

        pos = y[i,:]

        if np.isnan(pos[-1]):
            out.append(np.nan)
            continue
        fit = np.polyfit(t, pos, 1)
        out.append(fit[0])
    return np.array(out)


def measure(path):
    print path

    meas = dict()
    meas_names = ["speed","bound_cl","bound_mt"]

    for i, name in enumerate(meas_names):
        meas.update({name: i})
    dist = manual_parsing(os.path.join(path, "results.txt"), " ", float)
    couples = manual_parsing(os.path.join(path, "results_couples.txt"), " ", float)

    measurements = np.zeros([dist.shape[0],len(meas_names)])

    p = read_csv(os.path.join(path, "parameters.txt"), sep='|')
    t_max = p["time_step"] * p["nb_steps"]
    t = np.linspace(0, t_max[0], p["nb_frames"][0] + 1)
    measurements[:,meas["speed"]] = calc_speed(t[20:],dist[:,20:])
    measurements[:,meas["bound_cl"]] = couples[:,0]
    measurements[:, meas["bound_mt"]] = couples[:,1]
    np.savetxt(os.path.join(path,"measurements.txt"), measurements, delimiter=' ', header=" ".join(meas_names),comments="")



def main(args):
    paths = []
    PARENT = sys.path[0] + '/'
    for arg in args:
        if os.path.isdir(arg):
            paths.append(PARENT+arg)
        else:
            sys.stderr.write("  Error: unexpected argument `%s'\n" % arg)
            sys.exit()
    if not paths:
        sys.stderr.write("  Error: you must specify directories\n")
        sys.exit()

    for p in paths:
        measure(p)
# ------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == 'help':
        print(__doc__)
    else:
        main(sys.argv[1:])