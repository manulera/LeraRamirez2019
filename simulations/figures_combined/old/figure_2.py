from plot_settings import *

def theory_unbinding(ratio,v0,fs,ku,stiffness):
    v1 = ku*fs/stiffness
    return 1./(1.+ratio*v0/v1)

def theory_diffusive(rho_m,rho_c,v0,fs,D):
    v1 = D*fs/0.0042
    return 1./(1.+rho_c/rho_m/(1.-rho_c)*v0/v1)

def theory_diffusive_nobinding(ratio,v0,fs,D):
    v1 = D*fs/0.0042
    return 1./(1.+ratio*v0/v1)


# Figure1D -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","diff_ocup","data_summary","runs3_2")
measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ')
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|')
var_color = "mt_number"

# Lattice unit
a = 0.008
# Length of the MT
L = 3.
v0 = parameters["mt_speed"]
fs = parameters["stall_force"]
rho_m  = measurements["bound_mts"]*a/L
rho_c  = measurements["bound_cls"]*a/L
D  = parameters["cl_diff"]
x_sims = rho_c/rho_m
y_sims = measurements["speed"]/v0/2.

for f, f_par in enumerate(np.unique(parameters[var_color])):
    rho_condition = a*f_par/L
    ind = parameters[var_color]==f_par
    ave_rho_m = np.mean(rho_m[ind])
    x_rho_c = np.linspace(0,1)
    x_theo = x_rho_c/ave_rho_m

    sc = plt.scatter(x_sims[ind], y_sims[ind],label="Simulations: $\\rho_m$ = "+str(rho_condition),alpha=0.6,marker=marker.next())
    theory = theory_diffusive(ave_rho_m,x_rho_c, v0[0], fs[0], D[0])
    plt.plot(x_theo, theory, label="Theory: $\\rho_m$ = "+str(rho_condition))

ratio = np.linspace(0,5.)
theory = theory_diffusive_nobinding(ratio,v0[0],fs[0],D[0])
plt.plot(ratio,theory,c="black",label="Theory: no occupancy")

plt.legend()
plt.ylabel("$v_f/2v_0$")
plt.xlabel("$\\rho_c/\\rho_m$")
plt.ylim([0,1])
plt.xlim([0,5])

plt.savefig(os.path.join("fig2","fig2D.pdf"), transparent=True)


# Figure1F -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","diff_ocup","data_summary","runs5_3")
measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ')
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|')
var_color = "cl_diff"

# Lattice unit
a = 0.008
# Length of the MT
L = 3.
v0 = parameters["mt_speed"]
fs = parameters["stall_force"]
rho_m  = measurements["bound_mts"]*a/L
rho_c  = measurements["bound_cls"]*a/L
D  = parameters["cl_diff"]
x_sims = rho_m
y_sims = measurements["speed"]/v0/2.

for f, f_par in enumerate(np.unique(parameters[var_color])):

    ind = parameters[var_color]==f_par

    x_theo = np.linspace(0,1)
    v1 = f_par * fs[0] / 0.0042
    ratio = v0[0]/v1
    label = "$\gamma_c/\gamma_m$=%.2f"%ratio

    sc = plt.scatter(x_sims[ind], y_sims[ind],label="Simulations, "+label,alpha=0.6,marker=marker.next())


    theory = theory_diffusive(x_theo,x_theo, v0[0], fs[0], f_par)
    plt.plot(x_theo, theory, label="Theory, "+label)


# theory = theory_diffusive_nobinding(1,v0[0],fs[0],x_theo)
# plt.plot(x_theo/0.0042*fs[0]/v0[0],theory,c="black",label="Theory: no occupancy")

# theory = theory_diffusive_nobinding(2,v0[0],fs[0],x_theo)
# plt.plot(x_theo,theory,c="black",label="Theory: no occupancy")

plt.legend(fontsize=10)
plt.ylabel("$v_f/2v_0$")
plt.xlabel("$\\rho$")
plt.ylim([0,1])
plt.xlim([0,1])
plt.savefig(os.path.join("fig2","fig2F.pdf"), transparent=True)

plt.show()