from plot_settings import *
from theory_equations import *
import os
from pandas import read_csv

# Figure1A -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","Fig2A_bivalent_unbinding","data_summary","runs")
measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)
var_color = "mt_number"
x_sims = measurements["bound_cl"]/measurements["bound_mt"]
y_sims = measurements["speed"]/parameters["mt_speed"]/2.
v0 = parameters["mt_speed"]
fs = parameters["stall_force"]
ku = parameters["cl_unbinding"]
stiff = parameters["cl_stiffness"]


for f, f_par in enumerate(np.unique(parameters[var_color])):
    ind = parameters[var_color]==f_par
    sc = plt.scatter(x_sims[ind], y_sims[ind],label="Simulations: "+str(int(f_par))+" motors",alpha=0.6,marker=marker.next())

ratio = np.linspace(0,5)

theory = theory_bivalent_unbinding(ratio,v0[0],fs[0],ku[0],stiff[0])
plt.plot(ratio,theory,c="black",label="Theory")
plt.legend(frameon=False)
plt.ylabel("$v_{fil}\\ /\\ v_0$")
plt.xlabel("$c/m$")
plt.ylim([0,1])
plt.xlim([0,5])
plt.savefig(os.path.join("fig2","fig2A.pdf"), transparent=True)


# Figure1D -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","Fig2BC_bivalent_crosslinkers","data_summary","runs1")
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

legend_cols=["$\\rho_m$ ","Theory","Simul. "]
legend_rows=list()
scatters = list()
plots = list()

for f, f_par in enumerate(np.unique(parameters[var_color])):
    rho_condition = a*f_par/L
    ind = parameters[var_color]==f_par
    ave_rho_m = np.mean(rho_m[ind])
    x_rho_c = np.linspace(0.01,1)
    x_theo = x_rho_c/ave_rho_m

    scatters.append(plt.scatter(x_sims[ind], y_sims[ind],alpha=0.6,marker=marker.next()))
    legend_rows.append(str(rho_condition))
    theory = list()
    for rho_c_i in x_rho_c:
        theory.append(theory_bivalent_diffusive_continuous(ave_rho_m,rho_c_i, v0[0], fs[0], D[0]))

    plots.append(plt.plot(x_theo, theory, label="Theory: $\\rho_m$ = "+str(rho_condition)))

ratio = np.linspace(0,5.)
theory = theory_bivalent_diffusive_continuous_nobinding(ratio,v0[0],fs[0],D[0])
plt.plot(ratio,theory,c="black",label="Theory: no occupancy")

table_legend([plots,scatters],legend_rows,legend_cols,title_space=3,title_prespace=9)

plt.ylabel("$v_{fil}\\ /\\ v_0$")
plt.xlabel("$\\rho_c/\\rho_m$")
plt.ylim([0,1])
plt.xlim([0,5])

plt.savefig(os.path.join("fig2","fig2B.pdf"), transparent=True)

# Figure1F -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","Fig2BC_bivalent_crosslinkers","data_summary","runs2")
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

legend_cols=["$\gamma_c/\gamma_m$","Theory","Simul. "]
legend_rows=list()
scatters = list()
plots = list()

for f, f_par in enumerate(np.unique(parameters[var_color])):

    ind = parameters[var_color]==f_par

    x_theo = np.linspace(0.01,1)
    v1 = f_par * fs[0] / 0.0042
    ratio = v0[0]/v1
    legend_rows.append("%.2f"%ratio)

    scatters.append(plt.scatter(x_sims[ind], y_sims[ind],alpha=0.6,marker=marker.next()))
    theory = list()
    for xx in x_theo:
        theory.append(theory_bivalent_diffusive_continuous(xx,xx, v0[0], fs[0], f_par))
    plots.append(plt.plot(x_theo, theory))

table_legend([plots,scatters],legend_rows,legend_cols,title_space=3,title_prespace=8,bbox_to_anchor=(0.65,0.27))

plt.legend(fontsize=10)
plt.ylabel("$v_{fil}\\ /\\ v_0$")
plt.xlabel("$\\rho$")
plt.ylim([0,1])
plt.xlim([0,1])
plt.savefig(os.path.join("fig2","fig2C.pdf"), transparent=True)

plt.show()