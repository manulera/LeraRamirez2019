
#!/usr/bin/env python

import sys, os
import numpy as np
from pandas import read_csv,DataFrame

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


def measure(path):
    print path

    params = read_csv(os.path.join(path, "parameters.txt"), sep='|')
    meas_names = ["force", "bound_mts"]
    measurements = np.zeros([params.shape[0], len(meas_names)])
    measurements = DataFrame(measurements, columns=meas_names)

    force = manual_parsing(os.path.join(path, "results.txt"), " ", float)

    couples = manual_parsing(os.path.join(path, "results_couples.txt"), " ", float)


    measurements["bound_mts"] = np.mean(couples[:,20:],axis=1)
    measurements["force"] = np.mean(force[:, 20:]/couples[:,20:], axis=1)
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