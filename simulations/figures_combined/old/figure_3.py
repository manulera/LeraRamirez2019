from plot_settings import *


kT = 0.0042
def theory_force_diff(ratio):
    return 1./(1+ratio)

def theory_speed_diff(L,visc,D,v0,m,fs):

    drag_tail = 0.0042/D
    drag_MT = 3 * L / (np.log(0.5 * L / 0.0125) + 0.312) * np.pi*visc
    return v0/(1 + drag_MT*(1./drag_tail+v0/fs)/(2*m))

def theory_speed_diff_cls(Dc,Dm,v0,fs,rho_c,rho_m):
    term1 = rho_c/2./rho_m/(1.-rho_c)
    term2 = kT/Dc*(Dm/kT+v0/fs)
    return 1./(1.+term1*term2)

def speed2(Dc,Dm,v0,fs,rho_c,rho_m,L):

    term1 = v0+Dm/rho_m/L*np.log(1-rho_c)
    return term1*theory_speed_diff_cls(Dc,Dm,v0,fs,rho_c,rho_m)/v0

    # # ent = 0.0042 / 0.008 * rho_c
    # drag_m = kT/Dm
    # term1 = rho_c / rho_m / (1. - rho_c)
    # term2 = kT / Dc * (Dm / kT + v0 / fs)
    # if rho_m==0:
    #     ent=0
    # else:
    #     ent = ent / m / drag_m
    # return -(v0-ent)/(1 + term1 * term2)



# Figure2A -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","diffusive_motor_force_new","data_summary","runs1_new")
measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)
var_color = "mt_speed"
D = parameters["cl_diff"]
v0 = parameters["mt_speed"]
fs = parameters["stall_force"]

x_sims =fs*D/v0/kT
y_sims = measurements["force"]/fs/2.

for f, f_par in enumerate(np.unique(parameters[var_color])):
    ind = parameters[var_color]==f_par
    this_v0 = np.array(v0[ind])[0]
    sc = plt.scatter(x_sims[ind], y_sims[ind],label="Simulations",alpha=0.6,marker=marker.next())
    # sort_plot(parameters[var_ax][ind], predictions["bound_mt"][ind],plt.plot,color=sc.get_facecolors()[0])

ratio = np.logspace(-3, 3)
theory = theory_force_diff(ratio)
plt.semilogx(ratio,theory,color="black",label="Theory")

ratio = np.logspace(-3, 3)
plt.semilogx(ratio,1./(ratio/fs[0])/fs[0],color="black",label="$\gamma_d/\gamma_m$",linestyle="--")

# plt.plot(ratio,theory,c="black",label="Theory")
plt.legend(loc='lower left',frameon=False)
plt.ylabel("$f/f_s$")
plt.xlabel("$\gamma_m/\gamma_d$")
plt.ylim([0,1])
plt.xlim([1e-3,1e3])
plt.savefig(os.path.join("fig3","fig3A.pdf"), transparent=True)

# Figure2B -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","diffusive_motor_sliding","data_summary","runs1_new")

measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)

var_color = "viscosity"

D = parameters["cl_diff"]
v0 = parameters["mt_speed"]
fs = parameters["stall_force"]
L = parameters["fiber_length"]
motors = measurements["bound_mts"]


x_sims =motors/L
y_sims = measurements["speed"]

for f, f_par in enumerate(np.unique(parameters[var_color])):
    ind = parameters[var_color]==f_par

    sc = plt.scatter(x_sims[ind], y_sims[ind],label="Simulations: $\\xi$ "+str(f_par)+" $Pa.s$",alpha=0.6,marker=marker.next())
    dens_theo = np.linspace(0,15,100)
    theory = theory_speed_diff(L[0], f_par, D[0], v0[0], dens_theo*L[0], fs[0])

    plt.plot(dens_theo,theory,label="Theory: $\\xi$ "+str(f_par)+" $Pa.s$")
    # sort_plot(parameters[var_ax][ind], predictions["bound_mt"][ind],plt.plot,color=sc.get_facecolors()[0])

plt.xlabel("Motor density (motors/$\mu m$)")
plt.ylabel("Sliding speed ($\mu m/s$)")
plt.ylim([0,0.3])
plt.yticks([0,0.1,0.2,0.3])
plt.xlim([0,15])
plt.legend(fontsize=12,ncol=2)
plt.savefig(os.path.join("fig3","fig3B.pdf"), transparent=True)

# Figure2D -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","diffusive_motor_sliding_crosslinkers","data_summary","runs1_lat")

measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)

var_color = "cl_number"

Dc = parameters["cl_diff"]
Dm = parameters["mt_diff"]

v0 = parameters["mt_speed"]
fs = parameters["stall_force"]
L = parameters["fiber_length"]

rho_m = measurements["bound_mts"]*0.008/L
rho_c = parameters["cl_number"]*0.008/L


x_sims =rho_m
y_sims = measurements["speed"]/v0[0]

for f, f_par in enumerate(np.unique(parameters[var_color])):
    ind = parameters[var_color]==f_par
    ave_rho_c = np.array(rho_c[ind])[0]
    sc = plt.scatter(x_sims[ind], y_sims[ind],label="$\\rho_c$: %.2f" % ave_rho_c,alpha=0.6,marker=marker.next())


    rho_m_theo = np.linspace(0,0.5)
    theory = speed2(Dc[0],Dm[0],v0[0],fs[0],ave_rho_c,rho_m_theo,L[0])
    plt.plot(rho_m_theo,theory)

    # sort_plot(parameters[var_ax][ind], predictions["bound_mt"][ind],plt.plot,color=sc.get_facecolors()[0])

plt.xlabel("$\\rho_m$")
plt.ylabel("$v_0/v_f$")
plt.ylim([0,1])
plt.xlim([0,0.4])
# plt.yticks([0,0.1,0.2,0.3])
# plt.xlim([0,15])
plt.legend()
plt.savefig(os.path.join("fig3","fig3D.pdf"), transparent=True)

# Figure2E -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","diffusive_motor_sliding_crosslinkers","data_summary","runs_3_cls3")

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
extra_ax = plt.axes([.35, .5, .3, .3],facecolor='none')

plt.xticks([])
plt.yticks([])


for f, f_par in enumerate(np.unique(parameters[var_color])):
    plt.sca(main_axis)
    ind = parameters[var_color]==f_par
    ave_rho_m = np.array(rho_m[ind])[0]
    sc = plt.scatter(x_sims[ind], y_sims[ind],label="Simulations: $\\rho_m$: %.2f" % ave_rho_m,alpha=0.6,marker=marker.next())


    rho_c_theo = np.linspace(0,1)
    theory = speed2(Dc[0],Dm[0],v0[0],fs[0],rho_c_theo,ave_rho_m,L[0])
    theory2 = theory_speed_diff_cls(Dc[0],Dm[0],v0[0],fs[0],rho_c_theo,ave_rho_m)
    plt.plot(rho_c_theo,theory,label="Theory: $\\rho_m$: %.2f" % ave_rho_m)
    plt.axhline(y=0,linestyle='-.',c="grey",alpha=0.5)
    # sort_plot(parameters[var_ax][ind], predictions["bound_mt"][ind],plt.plot,color=sc.get_facecolors()[0])

    s = extra_ax.plot(rho_c_theo, theory,linewidth=1)

    extra_ax.plot(rho_c_theo, theory2,linestyle='--',c=s[0].get_color(),linewidth=1)

plt.xlabel("$\\rho_c$")
plt.ylabel("$v_0/v_f$")
# plt.ylim([0,1])
plt.xlim([0,1])
# extra_ax.axis['top'].set_visible(False)
# extra_ax.axis('off')
extra_ax.spines['right'].set_visible(False)
extra_ax.spines['top'].set_visible(False)
extra_ax.axhline(y=0,linestyle='-.',c="grey",alpha=0.5,linewidth=1)
extra_ax.set_xlim([0,1])
# extra_ax.axis['right'].set_visible(False)
# plt.yticks([0,0.1,0.2,0.3])
# plt.xlim([0,15])
# this is an inset axes over the main axes

plt.legend(fontsize=11)
plt.savefig(os.path.join("fig3","fig3E.pdf"), transparent=True)

plt.show()