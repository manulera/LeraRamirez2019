from plot_settings import *
from theory_equations import *
import os
from pandas import read_csv
import numpy as np


lansky_data = np.genfromtxt("data_lanksy.csv", dtype=float, delimiter=",", skip_header=6)

source_dir = ["runs1","runs2","runs3","runs4"]

for i,s in enumerate(source_dir):

    marker = new_fig()

    source = os.path.join("..","Fig4BC_entropic_expansion","data_summary",s)

    measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
    parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)

    x_sims =parameters["overlap_length"]
    y_sims = measurements["speed"]
    D = parameters["cl_diff"][0]

    occup = parameters["occup"][0]
    visc = parameters["viscosity"][0]
    stiffness = parameters["stiffness"][0]

    L_mt = 10.
    L = np.linspace(0,15,100)

    a = 0.008

    v= entropic_expansion_numerical(L,occup*L/a,D,visc,L_mt,stiffness)

    plt.scatter(x_sims,y_sims*1000,alpha=0.6,label="Simulations",marker=marker.next())
    plt.scatter(lansky_data[:, 0], lansky_data[:, 1],alpha=0.6 , label="Experimental Data",marker=marker.next())
    plt.plot(L,1000*v,label="Numer.",c="black",linestyle=":")
    v = entropic_expansion_continuous(L, occup * L / a, D, visc, L_mt, stiffness)
    plt.plot(L, 1000 * v, label="Analyt.", c="black")

    plt.ylim([-20,50])
    plt.xlim([0,15])
    plt.yticks(range(-10,51,10),[""]+map(str,range(0,50,10)))
    plt.axhline(y=0,linestyle='-.',c="black")
    plt.xlabel("Length ($\mu m$)")
    plt.ylabel("Expansion speed ($nm/s$)")
    plt.legend(frameon=False)
    plt.title("$D_1$: %.3f" % D + " $\mu m^2/s$, "+"$\kappa$: %u" % stiffness + " $pN/\mu m$, $\\rho_c$: %.2f" % occup,fontsize=18)
    plt.savefig("fig4/Fig4_%u.pdf" % stiffness)
plt.show()