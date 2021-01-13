
#!/usr/bin/env python

from equations import *

import sys, os
import numpy as np

from pandas import read_csv,DataFrame


def predict(path):
    print path


    params = read_csv(os.path.join(path, "parameters.txt"), sep='|',dtype=float)
    meas = read_csv(os.path.join(path, "measurements.txt"), sep=' ',dtype=float)

    pred_names = ["speed","bound_mt","bound_cl"]


    predictions = np.zeros([params.shape[0],len(pred_names)])
    predictions = DataFrame(predictions,columns=pred_names)

    v0 = params["mt_speed"]
    fs = params["stall_force"]
    kon_m = params["mt_binding"]
    koff_m = params["mt_unbinding"]
    kon_c = params["cl_binding"]
    koff_c = params["cl_unbinding"]
    mt = params["mt_number"]
    cl = params["cl_number"]
    stiffness = params["cl_stiffness"]
    predictions["bound_mt"] = binding_pred(kon_m,koff_m,mt,6.,v0)
    predictions["bound_cl"] = binding_pred(kon_c, koff_c, cl,6.,v0)
    predictions["speed"] = speed(v0,predictions["bound_cl"],predictions["bound_mt"],koff_c,stiffness,fs)
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