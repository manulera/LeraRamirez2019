from sympy import *

FF,AF,AA,kon,koff,v,L,T = var("FF AF AA kon koff v L T")

eq1 = -FF*kon*2 + AF*koff + AA*v/L*2
eq2 = FF*kon*2 - AF*(koff+kon) + AA*2*koff
eq3 = -AA*(2*koff+v/L*2) + AF*kon
eq4 = AA + AF + FF - T

solu = nonlinsolve([eq1,eq2,eq3,eq4],[AA,FF,AF,v])

print "a"