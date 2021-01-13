
from plot_settings import *
import numpy as np

kT = 0.0042
a = 0.008



def hop(f,stiffness):

    dG = f*a/kT - a*a/kT/2*stiffness
    if mode=="Lansky":
        return np.exp(dG/2.)
    elif mode=="Wang":
        return dG/(1-np.exp(-dG))





def drag(D,c,F,stiffness):
    f_m = F / c
    v_f = 2 * D / a * (hop(f_m, stiffness) - hop(-f_m, stiffness))
    return F/v_f



new_fig()
D = 0.1

c = np.linspace(0,400.)
force = 400.
for stiffness in [50.,100.,200.,300.]:

    mode = "Lansky"

    s = plt.plot(c/force, drag(D,c,force,stiffness), label="Lansky %u pN/$\mu$m" % stiffness, linewidth=3)

    mode = "Wang"

    plt.plot(c/force, drag(D, c, force, stiffness), label="Wang %u pN/$\mu$m" % stiffness, linewidth=3,color=s[0].get_color(),linestyle="--")
    # Method Lansky

drag_old = kT/D/2
plt.plot(c/force,drag_old*c,c="black",label="Continuous")
plt.xlabel("$c\\ /\\ F$ (pN$^{-1}$)")
plt.ylabel("Effective drag (pNs$/\mu m$)")

plt.legend(loc="best")
plt.xlim(xmin=0)
plt.ylim(ymin=0)
plt.savefig("fig4/fig4F.pdf")

plt.show()