import numpy as np
import matplotlib.pyplot as plt
from plot_settings import *
kT = 0.0042
mode = "Lansky"
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


new_fig()
plt.scatter(data[:,0],data[:,1],alpha=0.3,label="Lansky2015 data points")
L = np.linspace(0,40,100)
rho_c = 0.6
c = rho_c*L/a
D = 0.085
visc = 0.
L_mt = 0.
stiffness = 200.
plt.plot(L,1000*speed_entropic(L,c,D,visc,L_mt,stiffness),c="black",label="model")
plt.ylim([-10,50])
plt.ylabel("Expansion speed ($nm m/s$)")
plt.xlabel("Overlap length ($nm m/s$)")
plt.ylim(ymax=30)
plt.axhline(y=0,linestyle="-.",c="grey")
plt.tight_layout()
plt.legend()
plt.savefig("Lansky_fig.pdf")
plt.show()