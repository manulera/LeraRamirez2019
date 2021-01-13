import numpy as np
from Gilles import ReAct

# Initial conditions


# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)




def force(v0,D,fs,rho):
    f = 1./(D/0.0042/v0*rho+1./fs)
    speed = v0*(1-f/fs)
    return f,speed

def binding_pred(kon,koff,total,dist,L,S,v0,v):
    out = list()
    for i in range(total.shape[0]):
        out.append(solve_system(kon,koff,total[i],dist,L,S,v0[i],v[i]))
    return np.array(out)


def solve_system(kon,koff,total,dist,L,S,v0,v):


    L2 = L
    L1 = (S - L) / 2.
    L3 = L1

    user_input = ['FF', total,
                  'AF1', 0,
                  'AF2', 0,
                  'AF3', 0,
                  'AA', 0,
                  ]

    reactions = (

        (1, 'FF'), (1, 'AF1'), L1/S*kon,
        (1, 'AF1'), (1,'FF'),  koff,

        (1, 'FF'), (1, 'AF2'), L2/S*kon*2,
        (1, 'AF2'), (1, 'FF'), koff,

        (1, 'FF'), (1, 'AF3'), L3/S*kon,
        (1, 'AF3'), (1, 'FF'), koff,

        (1, 'AF2'), (1, 'AA'), kon,
        (1, 'AA'), (1, 'AF2'), 2*koff,

        (1, 'AA'), (1, 'AF1'), v/L2,

        (1, 'AF2'), (1, 'AA'), v0 / L2,

        (1, 'AF3'), (1, 'AF2'), v0 / L3 / 2,

    )

    solu = ReAct(user_input, reactions, 0, mode=3)[0]

    return solu[-1]