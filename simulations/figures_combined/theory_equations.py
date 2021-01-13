# This file contains the theoretical equations presented in Lera-Ramirez and Nedelec 2019.
# They are used for making the plots in the figures

from scipy.optimize import fsolve
import numpy as np

# Analytical solutions ---------------------------------------------------------

# Figure 2 000000000000000000000000000000000000000000000000


def theory_bivalent_unbinding(ratio,v0,fs,ku,stiffness):
    """
    The prediction for steady state speed in Fig. 1A (Crosslinker that unbind and dont diffuse)
    :param ratio: ratio of motors to crosslinkers
    :param v0: unloaded speed of the motor
    :param fs: stall force of the motor
    :param ku: unbinding rate of the crosslinkers
    :param stiffness: stiffness of the crosslinkers
    """
    v1 = ku*fs/stiffness
    return 1./(1.+ratio*v0/v1)


def theory_bivalent_diffusive_continuous(rho_m,rho_c,v0,fs,D,stiffness=None):
    """
    The prediction for steady state speed in Fig. 1B (Crosslinkers that diffuse and bivalent motors)
    :param rho_m: occupancy of motors (motors*lattice site length / overlap length)
    :param rho_c: occupancy of crosslinkers
    :param D: diffusion rate of crosslinker heads
    """
    v1 = D*fs/0.0042
    return 1./(1.+rho_c/rho_m/(1.-rho_c)*v0/v1)


def theory_bivalent_diffusive_continuous_nobinding(ratio,v0,fs,D,stiffness=None):
    """
    The prediction for steady state speed in Fig. 1B if occupancy was irrelevant.
    """
    v1 = D*fs/0.0042
    return 1./(1.+ratio*v0/v1)

# Figure 3 000000000000000000000000000000000000000000000000


def theory_force_diffussivemot_continuous(D,v0,fs,stiffness=None,kT=0.0042):
    """
    Returns the force exerted by a diffusive motor
    :param ratio: fs*D/v0/kT
    """
    ratio = D * fs / v0 / kT
    return 1./(1+ratio)


def theory_speed_diffussivemot_continuous(L,visc,D,v0,m,fs,stiffness=None):
    """
    Returns the steady state speed caused by diffusive motors
    :param L: Length of the microtubule
    :param visc: Viscosity of the medium
    :param D: Diffusion rate of the diffusive tail of the motor
    :param m: Number of motors
    """
    drag_tail = 0.0042/D
    drag_MT = 3 * L / (np.log(0.5 * L / 0.0125) + 0.312) * np.pi*visc
    return v0/(1 + drag_MT*(1./drag_tail+v0/fs)/(2*m))


def theory_speed_diffussivemot_diffusivecls_continuous(Dc,Dm,v0,fs,rho_c,rho_m,L,stiffness=None,kT=0.0042):
    """
    Returns the steady state sliding speed for diffusive motors and diffusive crosslinkers
    :param Dc: Diffusive rate of the crosslinker
    :param Dm: Diffusive rate of the motor tail
    """

    # Terms independent of entropic forces
    term1 = rho_c/2./rho_m/(1.-rho_c)
    term2 = kT/Dc*(Dm/kT+v0/fs)

    # Entropic forces term
    term3 = v0 + Dm / rho_m / L * np.log(1 - rho_c)
    return term3/(1.+term1*term2)/v0

# Figure 4 000000000000000000000000000000000000000000000000


def final_length_entropic_continuous(Dm,rho_m,v0,Lc,L_max,fs,stiffness=None):
    """
    Steady state overlap length reached for a given number of motors and crosslinkers
    :param Dm:
    :param rho_m:
    :param v0:
    :param Lc: Length occupied by all crosslinkers (c*a)
    :param L_max: Maximum overlap possible (Limited by the length of microtubules)
    :return:
    """
    def eq2solve(L, Dm, rho_m, v0, Lc):
        return rho_m * L * v0 / Dm + np.log(1 - Lc / L)
    out = list()
    for i in range(Lc.shape[0]):
        res = fsolve(eq2solve, Lc[i]*1.5, args=(Dm, rho_m, v0, Lc[i]))[0]
        if res > L_max:
            out.append(L_max)
        else:
            out.append(res)

    return np.array(out)


def entropic_expansion_continuous(L,c,D,visc,L_mt,stiffness=None,a=0.008,kT=0.0042):
    """
    Speed of expansion caused by entropic forces when only crosslinkers are present.
    :param L: Overlap length
    :param L_mt: Length of the microtubule (negligible, but for high viscosities would matter)
    """
    rho_c = c * a / L

    entropic_force = -kT / a * np.log(1 - c * a / L)

    # One can ignore this for small viscosities
    drag_MT = 3 * L_mt / (np.log(0.5 * L_mt / 0.0125) + 0.312) * np.pi * visc

    drag_cls = c * kT / D / (1. - rho_c) / 2.
    drag = drag_cls + drag_MT

    return entropic_force / drag


# Numerical solutions ---------------------------------------------------------


# General functions 000000000000000000000000000000000000000000000000

def hop(f,stiffness,a=0.008,kT=0.0042,mode="Wang"):
    """
    Calculates the force dependent term of the rates, using the Wang (dG/(1-np.exp(-dG))) or Lansky (exp(dG/2.))
    formulation

    :param f: Force applied on the crosslinker head
    :param stiffness:
    :param a: lattice site (default the lattice site of the microtubule)
    :param kT:
    :param mode: Formulation of the dependency of the rates on the force, as proposed by Lansky et al. 2015, or by our
    paper, based on Wang et al. 2003
    """

    dG = f*a/kT - a*a/kT/2*stiffness
    if mode=="Lansky":
        return np.exp(dG/2.)
    elif mode=="Wang":
        out = dG/(1-np.exp(-dG))
        ind = np.equal(dG, 0)
        out[ind]=1
        return out


def speed_cl(D,f,stiffness,a=0.008):
    """
    Average speed of the crosslinker under a constant force "f"
    """
    return D / a * (hop(f, stiffness) - hop(-f, stiffness))

# Figure 2 000000000000000000000000000000000000000000000000


def theory_bivalent_diffusive_numerical(rho_m,rho_c,v0,fs,D,stiffness,a=0.008):
    """
    Same as 'theory_diffusive_continuous', but solved numerically
    """
    def eq2solve(v_f, rho_m, rho_c, v0, fs, D, stiffness):
        gamma_m = fs / v0
        f_m = (v0 - v_f) * gamma_m
        f_c = f_m * rho_m / rho_c
        return -v_f + speed_cl(D, f_c, stiffness) * (1 - rho_c)

    return fsolve(eq2solve,0,args=(rho_m,rho_c,v0,fs,D,stiffness))/v0


def theory_bivalent_diffusive_numerical_nobinding(ratio,v0,fs,D,stiffness,a=0.008):
    """
        Same as 'theory_diffusive_continuous_nobinding', but solved numerically
    """
    def eq2solve(v_f, ratio, v0, fs, D, stiffness):
        gamma_m = fs / v0
        f_m = (v0 - v_f) * gamma_m
        f_c = f_m / ratio
        return -v_f + speed_cl(D, f_c, stiffness)

    return fsolve(eq2solve,0,args=(ratio,v0,fs,D,stiffness))/v0

# Figure 3 000000000000000000000000000000000000000000000000


def theory_force_diffussivemot_numerical(D,v0,fs,stiffness):
    """
    Same as 'theory_force_diffussivemot_continuous', but numerical
    """

    def eq2solve(f_m, D, v0, fs, stiffness):
        return speed_cl(D, f_m, stiffness) - v0 * (1 - f_m / fs)

    return fsolve(eq2solve, v0*D/0.0042, args=(D, v0, fs, stiffness)) / fs


def theory_speed_diffussivemot_numerical(L,visc,D,v0,m,fs,stiffness):
    """
    Same as 'theory_speed_diffussivemot_continuous', but numerical solution
    """
    def eq2solve(v_f, L, visc, D, v0, m, fs, stiffness):
        drag_MT = 3 * L / (np.log(0.5 * L / 0.0125) + 0.312) * np.pi * visc
        v_m = v0 - drag_MT * v0 / fs * v_f / m
        f_m = drag_MT * v_f / m
        return -v_m + 2 * v_f + speed_cl(D, f_m, stiffness)

    return fsolve(eq2solve,0,args=(L,visc,D,v0,m,fs,stiffness))*2


def theory_speed_diffussivemot_diffusivecls_numerical(Dc,Dm,v0,fs,rho_c,rho_m,L,stiffness,kT=0.0042):
    """
    Same as 'theory_speed_diffussivemot_diffusivecls_continuous', but numerical
    """
    def eq2solve(f_c, Dc, Dm, v0, fs, rho_c, rho_m, L, stiffness):
        f_m = rho_c * f_c / rho_m - kT / rho_m / L * np.log(1 - rho_c)
        vm = v0 * (1 - f_m / fs)
        vt = speed_cl(Dm, f_m, stiffness)
        vc = speed_cl(Dc, f_c, stiffness) * (1 - rho_c)
        return 2 * vc - vm + vt

    f_c = fsolve(eq2solve, 0., args=(Dc, Dm, v0, fs, rho_c, rho_m, L, stiffness))
    return speed_cl(Dc, f_c, stiffness)*(1-rho_c)*2/v0


# Figure 4 000000000000000000000000000000000000000000000000


def final_length_entropic_numerical(Dm,rho_m,v0,Lc,L_max,fs,stiffness,kT=0.0042,a=0.008):
    """
    Same as 'final_length_entropic_continuous', numerically
    """

    def eq2solve_1(f_m, Dm, rho_m, v0, Lc):
        # One could potentially remove the fm/fs term since its always way smaller than 1.
        return v0*(1-f_m/fs) - speed_cl(Dm, f_m, stiffness)

    def eq2solve_2(L, Dm, rho_m, v0, Lc, f_m):
        return rho_m * L / a * f_m + kT / a * np.log(1 - Lc / L)

    out = list()

    for i in range(Lc.shape[0]):

        f_m = fsolve(eq2solve_1, 0.1, args=(Dm, rho_m, v0, Lc[i]))[0]
        res = fsolve(eq2solve_2, 3, args=(Dm, rho_m, v0, Lc[i], f_m))

        if res > L_max:
            out.append(L_max)
        else:
            out.append(res)

    return np.array(out)


def entropic_expansion_numerical(L,c,D,visc,L_mt,stiffness,a=0.008,kT=0.0042):
    """
    Same as 'entropic_expansion_continuous', but using the forward and backward rates.
    The drag term is not included, since it is negligible for the cases considered in the paper.
    """
    f_c = -kT / a * np.log(1 - c * a / L) / c
    return (1 - c * a / L) * speed_cl(D, f_c, stiffness) * 2
