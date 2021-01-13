import numpy as np
from Gilles import ReAct

# Initial conditions


# Reaction template ((stoch_1,reactant_1,stoch_2,reactant_2),(stoch_1,product_1,stoch_2,product_2),k)




def binding_pred(kon,koff,total,L,v0):
    out = list()
    for i in range(total.shape[0]):
        out.append(solve_system(kon,koff,total[i],L[i],v0[i]))
    return np.array(out)


def solve_system(kon,koff,total,L,v0):

    user_input = ['FF', total,
                  'AF', 0,
                  'AA', 0,
                  ]

    reactions = (

        (1, 'FF'), (1, 'AF'), kon*2,
        (1, 'AF'), (1, 'FF'), koff,

        (1, 'AF'), (1, 'AA'), kon,
        (1, 'AA'), (1, 'AF'), 2*koff,

        (1, 'AF'), (1, 'FF'), 2 * v0/L,
        (1, 'AA'), (1, 'FF'), 2 * v0 / L,
    )

    solu = ReAct(user_input, reactions, 0, mode=3)[0]

    return solu[-1]

def speed(L,visc,D,v0,m,fs):

    drag_tail = 0.0042/D
    drag_MT = 3 * L / (np.log(0.5 * L / 0.0125) + 0.312) * np.pi*visc
    return v0/(1 + drag_MT*(1./drag_tail+v0/fs)/(2*m))