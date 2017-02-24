from __future__ import division
from sympy import *
from sympy.plotting import plot

from sympy.functions import exp
from sympy.core.containers import Tuple

x,y,rs,theta = symbols('x y r \theta',real=True)
c = symbols('c')
r = sqrt(x**2+y**2)
rp = sqrt(cos(theta)**2+y**2)
phi = symbols(r'\phi',cls=Function)

print '---------------------------------------'
print '           u (i = 1, l = 1)'
print '---------------------------------------'
u = diff(diff(phi(r),x),x)+diff(diff(phi(r),y),y) - diff(diff(phi(r),x),x)
u = simplify(u)
pprint(expand(u))

print '---------------------------------------'
print '           u (i = 1, l = 2)'
print '---------------------------------------'
u = -diff(diff(phi(r),x),y)
#u = simplify(u).subs(x**2+y**2,rs**2)
u = simplify(u)
pprint(expand(u))

print '---------------------------------------'
print '           u (i = 2, l = 1)'
print '---------------------------------------'
u = -diff(diff(phi(r),y),x)
u = simplify(u)
pprint(expand(u))

print '---------------------------------------'
print '           u (i = 2, l = 2)'
print '---------------------------------------'
u = diff(diff(phi(r),x),x)+diff(diff(phi(r),y),y) - diff(diff(phi(r),y),y)
u = simplify(u)
pprint(expand(u))

print '---------------------------------------'
print '           p (l = 1) '
print '---------------------------------------'
p = diff(diff(diff(phi(r),x),x),x)+diff(diff(diff(phi(r),y),y),x)
p = simplify(p)
pprint(simplify(expand(p)))

print '---------------------------------------'
print '           p (l = 2) '
print '---------------------------------------'
p = diff(diff(diff(phi(r),x),x),y)+diff(diff(diff(phi(r),y),y),y)
p = simplify(p)
pprint(simplify(expand(p)))

print '---------------------------------------'
print '           dudx (i = 1, l = 1)'
print '---------------------------------------'
u = diff(diff(phi(r),x),x)+diff(diff(phi(r),y),y) - diff(diff(phi(r),x),x)
u = diff(u,x)
u = simplify(u)
pprint(simplify(expand(u)))

print '---------------------------------------'
print '           dudx (i = 1, l = 2)'
print '---------------------------------------'
u = -diff(diff(phi(r),x),y)
u = diff(u,x)
u = simplify(u)
pprint(simplify(expand(u)))

print '---------------------------------------'
print '           dudx (i = 2, l = 1)'
print '---------------------------------------'
u = -diff(diff(phi(r),y),x)
u = diff(u,x)
u = simplify(u)
pprint(simplify(expand(u)))

print '---------------------------------------'
print '           dudx (i = 2, l = 2)'
print '---------------------------------------'
u = diff(diff(phi(r),x),x)+diff(diff(phi(r),y),y) - diff(diff(phi(r),y),y)
u = diff(u,x)
u = simplify(u)
pprint(simplify(expand(u)))

print '---------------------------------------'
print '           dudy (i = 1, l = 1)'
print '---------------------------------------'
u = diff(diff(phi(r),x),x)+diff(diff(phi(r),y),y) - diff(diff(phi(r),x),x)
u = diff(u,x)
u = simplify(u)
pprint(simplify(expand(u)))

print '---------------------------------------'
print '           dudy (i = 1, l = 2)'
print '---------------------------------------'
u = -diff(diff(phi(r),x),y)
u = diff(u,x)
u = simplify(u)
pprint(simplify(expand(u)))

print '---------------------------------------'
print '           dudy (i = 2, l = 1)'
print '---------------------------------------'
u = -diff(diff(phi(r),y),x)
u = diff(u,x)
u = simplify(u)
pprint(simplify(expand(u)))

print '---------------------------------------'
print '           dudy (i = 2, l = 2)'
print '---------------------------------------'
u = diff(diff(phi(r),x),x)+diff(diff(phi(r),y),y) - diff(diff(phi(r),y),y)
u = diff(u,x)
u = simplify(u)
pprint(simplify(expand(u)))



print '---------------------------------------'
phii = (sqrt((r**2+c**2))*(4*r**2+48*r**2*c**2-61*c**4)/75 - c**3*log(c)*r**2 - (5*r**2-2*c**2)*c**3*log(c+sqrt(r**2+c**2))/5)/12
u1 = simplify(diff(diff(phii,y),y))
u2 = simplify(diff(diff(phii,y),x))
print '---------------------------------------'
v1 = simplify(diff(diff(phii,x),y))
v2 = simplify(diff(diff(phii,x),x))
print '---------------------------------------'
p1 = (diff(diff(diff(phii,x),x),x)+diff(diff(diff(phii,y),y),x))
p2 = (diff(diff(diff(phii,x),x),y)+diff(diff(diff(phii,y),y),y))
print '---------------------------------------'
particular = Matrix([[u1,u2],[v1,v2],[p1,p2]])
xval = 1.0
yval = 2.0
cval = 3.0
evaled = particular.subs([(x,xval),(y,yval),(c,cval)])
pprint(evaled)
pprint(evaled.rank())
