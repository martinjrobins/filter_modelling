from __future__ import division
from sympy import *
from sympy.plotting import plot

from sympy.functions import exp
from sympy.core.containers import Tuple

x,y,rs = symbols('x y r',real=True)
c = symbols('c')
r = sqrt(x**2+y**2)
phi = (1.0/(12.0))*((1.0/75.0)*sqrt(r**2+c**2)*(4*r**4+48*r**2*c**2-61*c**4)-c**3*log(c)*r**2-(1.0/5.0)*(5*r**2-2*c**2)*c**3*log(c+sqrt(r**2+c**2)))

print '---------------------------------------'
print '           u (i = 1, l = 1)'
print '---------------------------------------'
u = diff(diff(phi,x),x)+diff(diff(phi,y),y) + diff(diff(phi,x),x)
u = collect(u,sqrt(r**2+c**2))
pprint(u)


print '---------------------------------------'
print '           p '
print '---------------------------------------'
p = diff(diff(diff(phi,x),x),y)+diff(diff(diff(phi,y),y),y) + diff(diff(diff(phi,x),x),x)+diff(diff(diff(phi,y),y),x)
#p = simplify(p).subs(r,rs).subs(x**2+y**2,rs**2).subs(x**4+2*x**2*y**2+y**4,rs**4)
p = simplify(p)
pprint(p)
