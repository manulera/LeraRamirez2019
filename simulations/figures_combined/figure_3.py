from plot_settings import *
from theory_equations import *
import os
from pandas import read_csv

kT =0.0042

# Figure2A -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","Fig3A_diffusive_motor_force","data_summary","runs")
measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)
var_color = "couple_stiffness"
D = parameters["cl_diff"]
v0 = parameters["mt_speed"]
fs = parameters["stall_force"]

x_sims =D*fs[0]/v0[0]/kT
y_sims = measurements["force"]/fs/2.

for f, f_par in enumerate(np.unique(parameters[var_color])):
    ind = parameters[var_color]==f_par

    sc = plt.scatter(x_sims[ind], y_sims[ind],label="Simulations",alpha=0.6,marker=marker.next())

    D = np.logspace(-7, 0)
    theory = list()
    for D_i in D:
        theory.append(theory_force_diffussivemot_continuous(D_i,v0[0],fs[0],f_par))

    plt.semilogx(D * fs[0] / v0[0] / kT, theory, color="black", label="Theory")

ratio = np.logspace(-3, 3)

plt.semilogx(ratio,1./(ratio/fs[0])/fs[0],color="black",label="$\gamma_d/\gamma_m$",linestyle="--")

plt.legend(loc='lower left',frameon=False)
plt.ylabel("$f/f_s$")
plt.xlabel("$\gamma_m/\gamma_d$")
plt.ylim([0,1.05])
plt.xlim([1e-3,1e3])
plt.savefig(os.path.join("fig3","fig3A.pdf"), transparent=True)

# Figure2B -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","Fig3B_diffusive_motor","data_summary","runs")

measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)

var_color = "viscosity"

D = parameters["cl_diff"]
v0 = parameters["mt_speed"]
fs = parameters["stall_force"]
L = parameters["fiber_length"]
motors = measurements["bound_mts"]

legend_cols=["$\\xi$ ($Pa.s$)","Theory","Simul. "]
legend_rows=list()
scatters = list()
plots = list()


x_sims =motors/L
y_sims = measurements["speed"]

for f, f_par in enumerate(np.unique(parameters[var_color])):
    ind = parameters[var_color]==f_par

    legend_rows.append(str(f_par))
    scatters.append(plt.scatter(x_sims[ind], y_sims[ind],label="Simulations: $\\xi$ "+str(f_par),alpha=0.6,marker=marker.next()))
    dens_theo = np.linspace(0,15,100)
    theory = list()

    for d_i in dens_theo:
        theory.append(theory_speed_diffussivemot_continuous(L[0], f_par, D[0], v0[0], d_i * L[0], fs[0],100.))

    plots.append(plt.plot(dens_theo,theory))

plt.xlabel("Motor density (motors/$\mu m$)")
plt.ylabel("$2v_{fil}$ ($\mu m/s$)")
plt.ylim([0,0.3])
plt.yticks([0,0.1,0.2,0.3])
plt.xlim([0,15])
table_legend([plots,scatters],legend_rows,legend_cols,title_space=3,title_prespace=4)

plt.savefig(os.path.join("fig3","fig3B.pdf"), transparent=True)

# Figure2E -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","Fig3C_diffusive_motor_crosslinkers","data_summary","runs")

measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)

var_color = "mt_number"

Dc = parameters["cl_diff"]
Dm = parameters["mt_diff"]

v0 = parameters["mt_speed"]
fs = parameters["stall_force"]
L = parameters["fiber_length"]

rho_m = measurements["bound_mts"]*0.008/L
rho_c = parameters["cl_number"]*0.008/L


x_sims =rho_c
y_sims = measurements["speed"]/v0[0]

main_axis = plt.gca()

legend_cols=["$\\rho_m$ ","Theory","Simul."]
legend_rows=list()
scatters = list()
plots = list()

for f, f_par in enumerate(np.unique(parameters[var_color])):
    plt.sca(main_axis)
    ind = parameters[var_color]==f_par
    ave_rho_m = np.array(rho_m[ind])[0]
    scatters.append(plt.scatter(x_sims[ind], y_sims[ind],label="Simulations: $\\rho_m$: %.2f" % ave_rho_m,alpha=0.6,marker=marker.next()))


    rho_c_theo = np.linspace(0,1)
    theory = list()
    for rho_c_i in rho_c_theo:
        theory.append(theory_speed_diffussivemot_diffusivecls_continuous(Dc[0],Dm[0],v0[0],fs[0],rho_c_i,ave_rho_m,L[0],100.))

    plots.append(plt.plot(rho_c_theo,theory,label="Theory (ana): $\\rho_m$: %.2f" % ave_rho_m))
    legend_rows.append("%.2f" % ave_rho_m)

    plt.axhline(y=0,linestyle='-.',c="grey",alpha=0.5)


plt.xlabel("$\\rho_c$")
plt.ylabel("$2v_{fil}\\ /\\ v_0$")

plt.xlim([0,1])
table_legend([plots,scatters],legend_rows,legend_cols,title_space=4,title_prespace=4)

plt.savefig(os.path.join("fig3","fig3C.pdf"), transparent=True)

plt.show()