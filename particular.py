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
phi_sol = (sqrt(r**2+c**2)*(4*r**4+48*r**2*c**2-61*c**4)/75-c**3*log(c)*r**2-(5*r**2-2*c**2)*c**3*log(c+sqrt(r**2+c**2))/5)/12
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

print '---------------------------------------'
print '           phi\'\'\'\''
print '---------------------------------------'
phi_sol_dash_dash_dash_dash = diff(diff(diff(diff(phi_sol,r),r),r),r)
pprint(simplify(phi_sol_dash_dash_dash_dash))

print '---------------------------------------'
print '           phi\'\'\'\'\''
print '---------------------------------------'
phi_sol_dash_dash_dash_dash_dash = diff(diff(diff(diff(diff(phi_sol,r),r),r),r),r)
pprint(simplify(phi_sol_dash_dash_dash_dash_dash))




print '---------------------------------------'
print '           phi\' / r-phi\'\''
print '---------------------------------------'
test = phi_sol_dash_div_r - phi_sol_dash_dash
pprint(simplify(test))

print '---------------------------------------'
print '           phi\' / r+phi\'\''
print '---------------------------------------'
test = phi_sol_dash_div_r + phi_sol_dash_dash
pprint(simplify(test))

print '---------------------------------------'
print '           phi\' / r^2-phi\'\'/r'
print '---------------------------------------'
test = phi_sol_dash_div_r/r - phi_sol_dash_dash/r
pprint(simplify(test))

print '---------------------------------------'
print '           phi\' / r^2+phi\'\'/r'
print '---------------------------------------'
test = phi_sol_dash_div_r/r + phi_sol_dash_dash/r
pprint(simplify(test))





