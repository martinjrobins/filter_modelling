from __future__ import division
from sympy import *
from sympy.plotting import plot

from sympy.functions import exp
from sympy.core.containers import Tuple

r = symbols('r')
c = symbols('c')
phi = Function(r'\phi')(r)


print '---------------------------------------'
print '           phi'
print '---------------------------------------'
phi_sol = (1.0/(12.0))*((1.0/75.0)*sqrt(r**2+c**2)*(4*r**4+48*r**2*c**2-61*c**4)-c**3*log(c)*r**2-(1.0/5.0)*(5*r**2-2*c**2)*c**3*log(c+sqrt(r**2+c**2)))
pprint(simplify(phi_sol))

print '---------------------------------------'
print '           phi\' / r'
print '---------------------------------------'
phi_sol_dash_div_r = diff(phi_sol,r)/r
pprint(simplify(phi_sol_dash_div_r))

print '---------------------------------------'
print '           phi\'\''
print '---------------------------------------'
phi_sol_dash_dash = diff(diff(phi_sol,r),r)
pprint(simplify(phi_sol_dash_dash))

print '---------------------------------------'
print '           phi\'\'\''
print '---------------------------------------'
phi_sol_dash_dash_dash = diff(diff(diff(phi_sol,r),r),r)
pprint(simplify(phi_sol_dash_dash_dash))



