

import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

kT = 0.0042
mode = "Wang"
a = 0.008

# Possible combinations stiffness and diffusion that satisfy the diffusion

st = [0.0, 6.122448979591836, 12.244897959183673, 18.36734693877551, 24.489795918367346, 30.612244897959183, 36.73469387755102, 42.857142857142854, 48.97959183673469, 55.10204081632653, 61.224489795918366, 67.3469387755102, 73.46938775510203, 79.59183673469387, 85.71428571428571, 91.83673469387755, 97.95918367346938, 104.08163265306122, 110.20408163265306, 116.3265306122449, 122.44897959183673, 128.57142857142856, 134.6938775510204, 140.81632653061223, 146.93877551020407, 153.0612244897959, 159.18367346938774, 165.30612244897958, 171.42857142857142, 177.55102040816325, 183.6734693877551, 189.79591836734693, 195.91836734693877, 202.0408163265306, 208.16326530612244, 214.28571428571428, 220.40816326530611, 226.53061224489795, 232.6530612244898, 238.77551020408163, 244.89795918367346, 251.0204081632653, 257.1428571428571, 263.265306122449, 269.3877551020408, 275.51020408163265, 281.63265306122446, 287.7551020408163, 293.87755102040813, 300.0]
D  = [0.011, 0.021571428571428575, 0.023081632653061228, 0.023081632653061228, 0.023081632653061228, 0.02459183673469388, 0.02459183673469388, 0.027612244897959187, 0.02459183673469388, 0.026102040816326534, 0.02459183673469388, 0.027612244897959187, 0.027612244897959187, 0.027612244897959187, 0.027612244897959187, 0.027612244897959187, 0.02912244897959184, 0.02912244897959184, 0.030632653061224493, 0.030632653061224493, 0.030632653061224493, 0.02912244897959184, 0.030632653061224493, 0.02912244897959184, 0.03214285714285715, 0.03214285714285715, 0.03214285714285715, 0.0336530612244898, 0.036673469387755106, 0.0336530612244898, 0.0336530612244898, 0.03516326530612245, 0.03516326530612245, 0.036673469387755106, 0.03969387755102041, 0.03969387755102041, 0.03818367346938776, 0.041204081632653065, 0.03969387755102041, 0.03818367346938776, 0.041204081632653065, 0.045734693877551025, 0.045734693877551025, 0.04422448979591838, 0.04422448979591838, 0.04422448979591838, 0.045734693877551025, 0.045734693877551025, 0.05026530612244899, 0.05177551020408164]


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



# plt.scatter(data[:,0],data[:,1],alpha=0.3)

visc = 0.
L_mt = 0.

rho_c = 0.6
r2_all = list()

plt.figure()
plt.scatter(data[:, 0],data[:,1],label="Lansky points")
for D_i, st_i in zip(D,st):
    # print D_i,st_i

    L_data = data[:, 0]
    c = rho_c * L_data / a
    prediction = 1000*speed_entropic(L_data,c,D_i,visc,L_mt,st_i)
    r2_all.append(r2_score(data[:,1],prediction))


    L = np.linspace(0,40,100)
    c = rho_c * L / a

    prediction = 1000 * speed_entropic(L, c, D_i, visc, L_mt, st_i)
    plt.plot(L, prediction)


plt.plot(D,r2_all)
plt.xlabel("Length (um)")
plt.ylabel("Expansion speed (um/s)")
plt.legend()
plt.show()

