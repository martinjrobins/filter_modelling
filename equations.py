from __future__ import division
from sympy import *
from sympy.plotting import plot

from sympy.functions import exp
from sympy.core.containers import Tuple

#potential = sqrt(gradientx u ^2  +  gradient y u ^2) +
#            sqrt(gradientx v ^2  +  gradient y v ^2)

x = symbols(r'x',real=True)
y = symbols(r'y',real=True)
u = Function(r'u',real=True)(x,y)
v = Function(r'v',real=True)(x,y)


print '---------------------------------------'
print '           P'
print '---------------------------------------'

P = simplify(sqrt(diff(u,x)**2 + diff(u,y)**2) + sqrt(diff(v,x)**2 + diff(v,y)**2))
pprint(P)

print '---------------------------------------'
print '           dPdx'
print '---------------------------------------'
dPdx = diff(P,x)
pprint(dPdx)

print '---------------------------------------'
print '           dPdy'
print '---------------------------------------'
dPdy = diff(P,y)
pprint(dPdy)
