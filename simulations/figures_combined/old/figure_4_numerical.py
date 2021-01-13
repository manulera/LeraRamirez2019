from plot_settings import *
from scipy.optimize import fsolve

kT = 0.0042
mode = "Wang"
a = 0.008

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

def eq2solve_1(f_m,Dm,rho_m,v0,Lc):
    return v0-speed_cl(Dm,f_m,100.)

def eq2solve_2(L,Dm,rho_m,v0,Lc,f_m):
    return rho_m*L/a*f_m + kT/a*np.log(1-Lc/L)

def final_length(Dm,rho_m,v0,Lc,L_max):
    out = list()

    for i in range(Lc.shape[0]):

        f_m = fsolve(eq2solve_1,0.1,args=(Dm,rho_m,v0,Lc[i]))[0]
        res = fsolve(eq2solve_2,3,args=(Dm,rho_m,v0,Lc[i],f_m))
        if res>L_max:
            out.append(L_max)
        else:
            out.append(res)

    return np.array(out)

def drag(D,c,L,a,L_mt,visc,do_occup):
    rho_c = c*a/L
    if not do_occup:
        rho_c = 0
    drag_MT = 3 * L_mt / (np.log(0.5 * L_mt / 0.0125) + 0.312) * np.pi * visc

    return c*kT/D/(1.-rho_c)/2.+drag_MT

def entropic(c,L,a):
    l = L/a
    return -kT / a * np.log((l + 1 - c) / (l + 1))

def speed_entropic(L,c,D,visc,L_mt,stiffness):
    f_c = -kT/a*np.log(1-c*a/L)/c

    return (1-c*a/L)*speed_cl(D,f_c,stiffness)*2

def theory_speed_diff_cls(Dc,Dm,v0,fs,rho_c,rho_m):
    term1 = rho_c/2./rho_m/(1.-rho_c)
    term2 = kT/Dc*(Dm/kT+v0/fs)
    return 1./(1.+term1*term2)

def speed2(Dc,Dm,v0,fs,rho_c,rho_m,L,stiffness):
    f_c = fsolve(eq2solve_5, 0., args=(Dc, Dm, v0, fs, rho_c, rho_m, L, stiffness))
    return speed_cl(Dc, f_c, stiffness)*(1-rho_c)*2/v0

def eq2solve_5(f_c,Dc,Dm,v0,fs,rho_c,rho_m,L,stiffness):
    f_m = rho_c*f_c/rho_m - kT/rho_m/L*np.log(1-rho_c)
    vm = v0*(1-f_m/fs)
    vt = speed_cl(Dm,f_m,stiffness)
    vc = speed_cl(Dc, f_c, stiffness)*(1-rho_c)
    return 2*vc -vm +vt


marker = new_fig()

source = os.path.join("..","bivalent_dynamics","data_summary","runs_ss")
measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)
var_color="cl_number"
x_sims =parameters["mt_number"]
y_sims = measurements["final_length"]
v0 = parameters["mt_speed"]
Lc = parameters["cl_number"]*0.008

for f, f_par in enumerate(np.unique(Lc)):
    ind = Lc==f_par
    nb_cls = int(f_par/0.008)
    sc = plt.scatter(x_sims[ind], y_sims[ind],label="Simulations: "+str(nb_cls)+ " crosslinkers",alpha=0.6,marker=marker.next())
    plt.axhline(f_par,color=sc.get_facecolors()[0],alpha=1,label="Lc: "+str(nb_cls)+ " crosslinkers")


# x = np.linspace(0,0.5)
# plt.plot(x,x,color="black",label="Theory")
plt.legend(frameon=False,fontsize=10,ncol=2)
plt.ylabel("Overlap length ($\mu m$)")
plt.xlabel("Motors")
plt.ylim([0,3.5])

plt.savefig(os.path.join("fig4","fig4A.pdf"), transparent=True)


# Figure4B -----------------------------------------------------------------------------------------------------------
marker = new_fig()

source = os.path.join("..","diffusive_motor_sliding_crosslinkers_dynamics","data_summary","runs1_several_stiffness")

measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)
var_color="mt_number"
x_sims =parameters["cl_number"]
y_sims = measurements["final_length"]
rho_m = measurements["rho_m"]
Dm = parameters["cl_diff"]
L_max = parameters["fiber_length"]
v0 = parameters["mt_speed"]
Lc = x_sims*0.008


for f, f_par in enumerate(np.unique(parameters[var_color])):
    ind = parameters[var_color]==f_par
    ave_rho_m = np.mean(rho_m[ind])
    x_theo = np.linspace(0,200)

    y_theo = final_length(Dm[0],ave_rho_m,v0[0],x_theo*0.008,L_max[0])

    sc = plt.scatter(x_sims[ind], y_sims[ind],label="Simulations: " + str(int(f_par)) +" motors",alpha=0.6,marker=marker.next())
    plt.plot(x_theo,y_theo,label="Theory: " + str(int(f_par)) +" motors")
plt.plot(x_theo,x_theo*0.008,label="$L_c$", linestyle='--',c='black')

# x = np.linspace(0,0.5)
# plt.plot(x,x,color="black",label="Theory")
# plt.legend(frameon=False)
plt.ylabel("Overlap length ($\mu m$)")
plt.xlabel("Crosslinkers")
plt.ylim(ymin=0)
plt.xlim([0,200])
plt.legend()
plt.savefig(os.path.join("fig4","fig4B.pdf"), transparent=True)

# Figure4C -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source1 = os.path.join("..","bivalent_growth","data_summary","runs1")
source2 = os.path.join("..","diffusive_motor_sliding_crosslinkers_growth","data_summary","runs1_less_motors_new")

dirs = [source1,source2]
labels= ["bivalent", "diffusive"]
for i,d in enumerate(dirs):

    measurements = read_csv(os.path.join(d,"measurements.txt"), sep=' ',dtype=float)
    parameters = read_csv(os.path.join(d, "parameters.txt"), sep='|',dtype=float)

    x_sims =parameters["growth_speed"]/parameters["mt_speed"]
    y_sims = measurements["speed"]/parameters["mt_speed"]/2.
    c_num = parameters["cl_number"][0]
    if d==source2:
        L = parameters["fiber_length"][0]
        Dm = parameters["mt_diff"][0]
        Dc = parameters["cl_diff"][0]
        v0 = parameters["mt_speed"][0]
        ave_rhom = np.mean(measurements["bound_mts"][x_sims>0.4])*0.008/L

        # max_speed = speed2(Dc, Dm, v0, 6., c_num*0.008/L, ave_rhom, L,100.)/2.

         # = 1./(1+c_num/ave_m*Dm/Dc/2.)/2.
        # plt.plot([max_speed,0.8],[max_speed,max_speed],linestyle="--",c="black",label="Diffusive max speed")

    sc = plt.scatter(x_sims, y_sims,label="Simulations: "+ labels[i],alpha=0.6,marker=marker.next())


x = np.linspace(0,0.8)
plt.plot(x,x,color="black",label="$v_{fil}=v_g$")
plt.legend(frameon=False)
plt.ylabel("$v_{fil}/v_0$")
plt.xlabel("$v_g/v_0$")
plt.ylim([0,0.8])
plt.xlim([0,0.8])
plt.savefig(os.path.join("fig4","fig4C.pdf"), transparent=True)

# Figure4D -----------------------------------------------------------------------------------------------------------
marker = new_fig()

source = os.path.join("..","entropic_expansion","data_summary","runs2_visc001")
# source = os.path.join("..","entropic_expansion","data_summary","runs2_mmorevisc_higher_occ")
measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)

x_sims =parameters["overlap_length"]
y_sims = measurements["speed"]
D = parameters["cl_diff"][0]
D = D
occup = parameters["occup"]
visc = parameters["viscosity"][0]

L_mt = 10.
L = np.linspace(0,10,200)

their_curve = 0.085*2/L

# plt.plot(L,their_curve*1000,linestyle='--',c="black",label="Lanksy et al. 2015: linear friction")

a = 0.008

for f,f_par in enumerate(np.unique(occup)):
    ind = occup==f_par

    v= speed_entropic(L,f_par*L/a,D,visc,L_mt,100.)
    plt.scatter(x_sims[ind],y_sims[ind],alpha=0.6,label="Simulations")
    plt.plot(L,v,label="Theory: occupancy",c="black")
    # v, dra, ent = speed_entropic(L, f_par * L / a, D, L_mt, visc, a,do_occup=False)
    # plt.plot(L, v, label="Theory: no occupancy", c="black",linestyle="--")


plt.ylim([-0.02,0.2])
plt.yticks([0,0.1,0.2])
plt.axhline(y=0,linestyle='-.',c="grey",alpha=0.5)
plt.xlabel("Length ($\mu m$)")
plt.ylabel("Expansion speed ($\mu m/s$)")
plt.legend()
plt.savefig(os.path.join("fig4","fig4D.pdf"), transparent=True)

