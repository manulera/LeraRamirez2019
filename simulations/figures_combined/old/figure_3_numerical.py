from plot_settings import *
from scipy.optimize import fsolve

mode = "Wang"
a = 0.008
kT = 0.0042



def hop(f,stiffness):

    dG = f*a/kT - a*a/kT/2*stiffness
    if mode=="Lansky":
        return np.exp(dG/2.)
    elif mode=="Wang":
        return dG/(1-np.exp(-dG))

def speed_cl(D,f,stiffness):
    return D / a * (hop(f, stiffness) - hop(-f, stiffness))


def eq2solve_1(v_f,L,visc,D,v0,m,fs,stiffness):
    drag_MT = 3 * L / (np.log(0.5 * L / 0.0125) + 0.312) * np.pi * visc
    v_m = v0- drag_MT*v0/fs*v_f/m
    f_m = drag_MT * v_f/m
    return -v_m + 2*v_f + speed_cl(D,f_m,stiffness)

def theory_speed_diff(L,visc,D,v0,m,fs,stiffness):
    return fsolve(eq2solve_1,0,args=(L,visc,D,v0,m,fs,stiffness))*2

def eq2solve_2(f_m,D,v0,fs,stiffness):
    return speed_cl(D,f_m,stiffness) + v0*(f_m/fs-1)

def theory_force_diff(D,v0,fs,stiffness):
    return fsolve(eq2solve_2, fs/2, args=(D, v0, fs, stiffness))/fs


def speed2(Dc,Dm,v0,fs,rho_c,rho_m,L,stiffness):
    f_c = fsolve(eq2solve_5, 0., args=(Dc, Dm, v0, fs, rho_c, rho_m, L, stiffness))
    return speed_cl(Dc, f_c, stiffness)*(1-rho_c)*2/v0


def eq2solve_5(f_c,Dc,Dm,v0,fs,rho_c,rho_m,L,stiffness):
    f_m = rho_c*f_c/rho_m - kT/rho_m/L*np.log(1-rho_c)
    vm = v0*(1-f_m/fs)
    vt = speed_cl(Dm,f_m,stiffness)
    vc = speed_cl(Dc, f_c, stiffness)*(1-rho_c)
    return 2*vc -vm +vt



# Figure2A -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source = os.path.join("..","diffusive_motor_force_new","data_summary","runs1_new")
measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)
var_color = "mt_speed"
D = parameters["cl_diff"]
v0 = parameters["mt_speed"]
fs = parameters["stall_force"]

x_sims =D*fs[0]/v0[0]/kT
y_sims = measurements["force"]/fs/2.

for f, f_par in enumerate(np.unique(parameters[var_color])):
    ind = parameters[var_color]==f_par
    this_v0 = np.array(v0[ind])[0]
    sc = plt.scatter(x_sims[ind], y_sims[ind],label="Simulations",alpha=0.6,marker=marker.next())
    # sort_plot(parameters[var_ax][ind], predictions["bound_mt"][ind],plt.plot,color=sc.get_facecolors()[0])

ratio = np.logspace(-3, 3)
theory = list()
D = np.logspace(-7,0)
for D_i in D:
    theory.append(theory_force_diff(D_i,v0[0],fs[0],100.))

plt.semilogx(D*fs[0]/v0[0]/kT,theory,color="black",label="Theory")

# plt.semilogx(D*fs[0]/v0[0]/kT,kT*v0[0]/D,color="black",linestyle="--",label="$\gamma_d/\gamma_m$")
ratio = np.logspace(-3, 3)
plt.semilogx(ratio,1./(ratio/fs[0])/fs[0],color="black",label="$\gamma_d/\gamma_m$",linestyle="--")

# plt.plot(ratio,theory,c="black",label="Theory")
plt.legend(loc='lower left',frameon=False)
plt.ylabel("$f/f_s$")
plt.xlabel("$D$")
plt.ylim([0,1])
plt.xlim([1e-3,1e3])
plt.savefig(os.path.join("fig3","fig3A.pdf"), transparent=True)
plt.show()

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
    theory = list()
    for d_i in dens_theo:
        theory.append(theory_speed_diff(L[0], f_par, D[0], v0[0], d_i * L[0], fs[0], 100.))

    plt.plot(dens_theo,theory,label="Theory: $\\xi$ "+str(f_par)+" $Pa.s$")
    # sort_plot(parameters[var_ax][ind], predictions["bound_mt"][ind],plt.plot,color=sc.get_facecolors()[0])

plt.xlabel("Motor density (motors/$\mu m$)")
plt.ylabel("$2v_{fil}$ ($\mu m/s$)")
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
    theory = list()
    for rho_m_i in rho_m_theo:
        theory.append(speed2(Dc[0], Dm[0], v0[0], fs[0], ave_rho_c, rho_m_i, L[0], 100.))

    plt.plot(rho_m_theo,theory)

    # sort_plot(parameters[var_ax][ind], predictions["bound_mt"][ind],plt.plot,color=sc.get_facecolors()[0])

plt.xlabel("$\\rho_m$")
plt.ylabel("$2v_{fil}\\ /\\ v_0$")
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
    theory = list()
    for rho_c_i in rho_c_theo:
        theory.append(speed2(Dc[0],Dm[0],v0[0],fs[0],rho_c_i,ave_rho_m,L[0],100.))
    # theory2 = theory_speed_diff_cls(Dc[0],Dm[0],v0[0],fs[0],rho_c_theo,ave_rho_m)
    plt.plot(rho_c_theo,theory,label="Theory: $\\rho_m$: %.2f" % ave_rho_m)
    plt.axhline(y=0,linestyle='-.',c="grey",alpha=0.5)
    # sort_plot(parameters[var_ax][ind], predictions["bound_mt"][ind],plt.plot,color=sc.get_facecolors()[0])

    s = extra_ax.plot(rho_c_theo, theory,linewidth=1)

    # extra_ax.plot(rho_c_theo, theory2,linestyle='--',c=s[0].get_color(),linewidth=1)

plt.xlabel("$\\rho_c$")
plt.ylabel("$2v_{fil}\\ /\\ v_0$")
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