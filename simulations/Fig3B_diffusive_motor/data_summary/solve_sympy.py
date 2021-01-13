from sympy import solve
from sympy import var

Ldy, Ldz = var('Ldy Ldz')
g, x, y, z = var('g x y z')
xZ, yZ, zZ = var('xZ yZ zZ')
xdd, ydd, zdd = var('xdd ydd zdd')

E1 = z * xdd + (xZ - x) * (g + zdd)
E2 = z * ydd + (yZ - y) * (g + zdd) - Ldy
E3 = -y * xdd + x * ydd - zZ * (g + zdd) + Ldz

sols = solve([E1, E2, E3], [xdd, ydd, Ldy])

print "xdd = ", (sols[xdd]).factor()
print "ydd = ", (sols[ydd]).factor()
print "Ldy = ", (sols[Ldy]).factor()