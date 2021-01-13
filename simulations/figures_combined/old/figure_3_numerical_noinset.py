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
        out = dG/(1-np.exp(-dG))
        ind = np.equal(dG, 0)
        out[ind]=1
        return out

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

def speed2_noent(Dc,Dm,v0,fs,rho_c,rho_m,L,stiffness):
    f_c = fsolve(eq2solve_5_noent, 0., args=(Dc, Dm, v0, fs, rho_c, rho_m, L, stiffness))
    return speed_cl(Dc, f_c, stiffness)*(1-rho_c)*2/v0


def eq2solve_5_noent(f_c,Dc,Dm,v0,fs,rho_c,rho_m,L,stiffness):
    f_m = rho_c*f_c/rho_m
    vm = v0*(1-f_m/fs)
    vt = speed_cl(Dm,f_m,stiffness)
    vc = speed_cl(Dc, f_c, stiffness)*(1-rho_c)
    return 2*vc -vm +vt


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



for f, f_par in enumerate(np.unique(parameters[var_color])):
    plt.sca(main_axis)
    ind = parameters[var_color]==f_par
    ave_rho_m = np.array(rho_m[ind])[0]
    sc = plt.scatter(x_sims[ind], y_sims[ind],label="Simulations: $\\rho_m$: %.2f" % ave_rho_m,alpha=0.6,marker=marker.next())


    rho_c_theo = np.linspace(0,1)
    theory = list()
    for rho_c_i in rho_c_theo:
        theory.append(speed2(Dc[0],Dm[0],v0[0],fs[0],rho_c_i,ave_rho_m,L[0],100.))
    theory2 = list()
    for rho_c_i in rho_c_theo:
        theory2.append(speed2_noent(Dc[0],Dm[0],v0[0],fs[0],rho_c_i,ave_rho_m,L[0],100.))


    s = plt.plot(rho_c_theo,theory,label="Theory: $\\rho_m$: %.2f" % ave_rho_m)

    plt.plot(rho_c_theo, theory2, c=s[0].get_color(),linestyle="--",linewidth=1,alpha=0.8)

    plt.axhline(y=0,linestyle='-.',c="grey",alpha=0.5)




    # extra_ax.plot(rho_c_theo, theory2,linestyle='--',c=s[0].get_color(),linewidth=1)

plt.xlabel("$\\rho_c$")
plt.ylabel("$v_0/v_f$")
plt.xlim([0,1])
plt.ylim(ymax=1)

plt.legend(fontsize=11)
plt.savefig(os.path.join("fig3","fig3E_noinset.pdf"), transparent=True)

plt.show()