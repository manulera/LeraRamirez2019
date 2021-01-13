import numpy as np
from Gilles import ReAct

def speed(v0,c,m,ku,stiffness,fs):
    v1 = fs*ku/stiffness
    return 2*v0/(1.+(c-1)/m*v0/v1)



def binding_pred(kon,koff,total,L,v0):
    out = list()
    for i in range(total.shape[0]):
        out.append(solve_system(kon[i],koff[i],total[i],L,v0[i]))
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

        (1, 'AF'), (1, 'FF'), v0 / L,
        #(1, 'AA'), (1, 'FF'), 2 * v0 / L,
    )

    solu = ReAct(user_input, reactions, 0, mode=3)[0]

    return solu[-1]