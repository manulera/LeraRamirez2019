
#!/usr/bin/env python

from equations import *

import sys, os
import numpy as np

from pandas import read_csv,DataFrame


def predict(path):
    print path


    params = read_csv(os.path.join(path, "parameters.txt"), sep='|')
    meas = read_csv(os.path.join(path, "measurements.txt"), sep=' ')

    pred_names = ["force","force2","bound_mts","max_force"]


    predictions = np.zeros([params.shape[0],len(pred_names)])
    predictions = DataFrame(predictions,columns=pred_names)

    v0 = params["mt_speed"]
    D = params["cl_diff"]
    fs = params["stall_force"]


    total = params["mt_number"]
    kon = 1.
    koff = 1.
    dist = 0.05
    L = 6.

    rho=1
    pred_force,pred_speed  = force(v0, D,fs,rho)

    predictions["bound_mts"] = binding_pred(kon, koff, total, dist, L, 1, v0,pred_speed)

    predictions["force"] = pred_force
    omega = 0.3

    pred_force, pred_speed = force(v0, D, fs, omega)
    predictions["force2"] = pred_force

    predictions["max_force"] = fs

    np.savetxt(os.path.join(path,"predictions.txt"), np.array(predictions), delimiter=' ', header=" ".join(pred_names),comments="")



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
        predict(p)
# ------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] == 'help':
        print(__doc__)
    else:
        main(sys.argv[1:])