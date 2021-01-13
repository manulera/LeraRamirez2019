import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

kT = 0.0042
mode = "Wang"
a = 0.008

def hop(f,stiffness):

    dG = f*a/kT - a*a/kT/2*stiffness
    if mode=="Lansky":
        return np.exp(dG/2.)
    elif mode=="Wang":
        return dG/(1-np.exp(-dG))

def speed_cl(D,f,stiffness):
    return D / a * (hop(f, stiffness) - hop(-f, stiffness))

def speed_entropic(L,c,D,visc,L_mt,stiffness):
    f_c = -kT/a*np.log(1-c*a/L)/c

    return (1-c*a/L)*speed_cl(D,f_c,stiffness)*2


data = np.genfromtxt('lansky2015_Fig3B.csv',dtype=float,delimiter=",",skip_header=6)


plt.figure()
# plt.scatter(data[:,0],data[:,1],alpha=0.3)


D = 0.085
visc = 0.
L_mt = 0.
stiffness = 400.
rho_c_all = np.linspace(0.2,0.99)
r2_all = list()
for rho_c in rho_c_all:
    L_data = data[:, 0]
    c = rho_c * L_data / a
    prediction = 1000*speed_entropic(L_data,c,D,visc,L_mt,stiffness)
    r2_all.append(r2_score(data[:,1],prediction))


plt.plot(rho_c_all,r2_all)

plt.show()

