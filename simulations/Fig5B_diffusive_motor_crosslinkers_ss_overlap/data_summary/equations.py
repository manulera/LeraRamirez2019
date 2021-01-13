import numpy as np
from Gilles import ReAct
from scipy.optimize import fsolve
# Initial conditions


# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)




def binding_pred(kon,koff,total):
    out = list()
    for i in range(total.shape[0]):
        out.append(solve_system(kon,koff,total[i]))
    return np.array(out)


def solve_system(kon,koff,total):

    user_input = ['FF', total,
                  'AF', 0,
                  'AA', 0,
                  ]

    reactions = (

        (1, 'FF'), (1, 'AF'), kon*2,
        (1, 'AF'), (1, 'FF'), koff,

        (1, 'AF'), (1, 'AA'), kon,
        (1, 'AA'), (1, 'AF'), 2*koff,

    )

    solu = ReAct(user_input, reactions, 0, mode=3)[0]

    return solu[-1]

def speed(Dm,Dc,v0,m,c,L):
    rho_c = c*0.008/L
    rho_m = m * 0.008 / L
    return v0/(1 + Dm/Dc/rho_m*rho_c/(1-rho_c))

def eq2solve(L,Dm,rho_m,v0,Lc):
    return rho_m*L*v0/Dm + np.log(1-Lc/L)

def final_length(Dm,rho_m,v0,Lc,sugg=1.5):
    # if we pass m

    # return Lc/(1-np.exp(rho_m*0.008/Dm*v0))
    #
    out = list()
    for i in range(Dm.shape[0]):
        out.append(fsolve(eq2solve,Lc[i]*sugg,args=(Dm[i],rho_m[i],v0[i],Lc[i]))[0])
        print eq2solve(out[i], Dm[i], rho_m[i], v0[i], Lc[i])

    return np.array(out)