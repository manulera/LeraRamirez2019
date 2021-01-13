
import numpy as np

def speed(v0,D,fs,rho_c,rho_m,meas_speed):
    v1 = fs*D/0.0042
    x_show = rho_c / rho_m / (1. - rho_c)
    vm = 1. / (1. / v0 +x_show/v1)
    y_show = (v0-vm)*v1/vm/v0
    y_2 = (v0-meas_speed/2.)*v1/meas_speed/v0

    return vm*2,x_show,y_show,y_2
