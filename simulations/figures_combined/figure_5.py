from plot_settings import *
from theory_equations import *
import os
from pandas import read_csv

marker = new_fig()

source = os.path.join("..","Fig5A_bivalent_ss_overlap","data_summary","runs")
measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)
var_color="cl_number"
x_sims =parameters["mt_number"]
y_sims = measurements["final_length"]
v0 = parameters["mt_speed"]
Lc = parameters["cl_number"]*0.008

legend_cols=["       $c$   ","  $L_c$  ","   Simul."]
legend_rows=list()
scatters = list()
plots = list()


for f, f_par in enumerate(np.unique(Lc)):
    ind = Lc==f_par
    nb_cls = int(f_par/0.008)
    scatters.append(plt.scatter(x_sims[ind], y_sims[ind],label="Simulations: "+str(nb_cls)+ " crosslinkers",alpha=0.6,marker=marker.next()))
    plots.append(plt.axhline(f_par,color=scatters[-1].get_facecolors()[0],alpha=1,label="Lc: "+str(nb_cls)+ " crosslinkers"))
    legend_rows.append(str(nb_cls))


plt.ylabel("Overlap length ($\mu m$)")
plt.xlabel("Motors")
plt.ylim([0,3.5])
table_legend([plots,scatters],legend_rows,legend_cols,title_space=4,title_prespace=2)

plt.savefig(os.path.join("fig5","fig5A.pdf"), transparent=True)


# Figure4B -----------------------------------------------------------------------------------------------------------
marker = new_fig()

source = os.path.join("..","Fig5B_diffusive_motor_crosslinkers_ss_overlap","data_summary","runs")

measurements = read_csv(os.path.join(source,"measurements.txt"), sep=' ',dtype=float)
parameters = read_csv(os.path.join(source, "parameters.txt"), sep='|',dtype=float)
var_color="mt_number"
x_sims =parameters["cl_number"]
y_sims = measurements["final_length"]
rho_m = measurements["rho_m"]
Dm = parameters["cl_diff"]
L_max = parameters["fiber_length"]
v0 = parameters["mt_speed"]
fs = parameters["stall_force"]
Lc = x_sims*0.008

legend_cols=["$\\rho_m$","Theory","Simul."]


legend_rows=list()
scatters = list()
plots = list()
plots2 = list()

for f, f_par in enumerate(np.unique(parameters[var_color])):
    ind = parameters[var_color]==f_par
    ave_rho_m = np.mean(rho_m[ind])
    legend_rows.append("%.2f" % ave_rho_m)
    x_theo = np.linspace(0,200)

    y_theo = final_length_entropic_continuous(Dm[0],ave_rho_m,v0[0],x_theo*0.008,L_max[0],fs[0],100)
    scatters.append(plt.scatter(x_sims[ind], y_sims[ind],label="Simulations: " + str(int(f_par)) +" motors",alpha=0.6,marker=marker.next()))

    plots.append(plt.plot(x_theo,y_theo,label="Theory (ana): " + str(int(f_par)) +" motors"))

    y_theo = final_length_entropic_numerical(Dm[0], ave_rho_m, v0[0], x_theo * 0.008, L_max[0], fs[0], 100)

plt.plot(x_theo,x_theo*0.008,label="$L_c$", linestyle='--',c='black')


table_legend([plots,scatters],legend_rows,legend_cols,title_space=4,title_prespace=8,handletextpad=1)

plt.ylabel("Overlap length ($\mu m$)")
plt.xlabel("Crosslinkers")
plt.ylim(ymin=0)
plt.xlim([0,200])

plt.savefig(os.path.join("fig5","fig5B.pdf"), transparent=True)

# Figure4C -----------------------------------------------------------------------------------------------------------

marker = new_fig()

source1 = os.path.join("..","Fig5C_bivalent_growth","data_summary","runs")
source2 = os.path.join("..","Fig5C_diffusive_motor_crosslinkers_growth","data_summary","runs")

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


    sc = plt.scatter(x_sims, y_sims,label="Simulations: "+ labels[i],alpha=0.6,marker=marker.next())


x = np.linspace(0,0.8)
plt.plot(x,x,color="black",label="$v_{fil}=v_g$")
plt.legend(frameon=False)
plt.ylabel("$v_{fil}/v_0$")
plt.xlabel("$v_g/v_0$")
plt.ylim([0,0.8])
plt.xlim([0,0.8])
plt.savefig(os.path.join("fig5","fig5C.pdf"), transparent=True)

plt.show()
