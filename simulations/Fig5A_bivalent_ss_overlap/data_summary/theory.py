
#!/usr/bin/env python

from equations import *

import sys, os
import numpy as np

from pandas import read_csv,DataFrame


def predict(path):
    print path


    params = read_csv(os.path.join(path, "parameters.txt"), sep='|')
    meas = read_csv(os.path.join(path, "measurements.txt"), sep=' ')

    pred_names = ["speed","x_show","y_show","data_scaled"]


    predictions = np.zeros([params.shape[0],len(pred_names)])
    predictions = DataFrame(predictions,columns=pred_names)

    v0 = params["mt_speed"]
    D = params["cl_diff"]
    fs = params["stall_force"]
    rho_c = meas["bound_cls"]*0.008/3.
    rho_m = params["mt_number"]*0.008/3.
    predictions["speed"],predictions["x_show"],predictions["y_show"],predictions["data_scaled"] = speed(v0, D, fs, rho_c, rho_m,meas["speed"])

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