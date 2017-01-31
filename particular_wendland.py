from __future__ import division
from sympy import *
from sympy.plotting import plot

from sympy.functions import exp
from sympy.core.containers import Tuple

x,y,rs = symbols('x y r',real=True)
c = symbols('c')
r = sqrt(x**2+y**2)
phi = symbols(r'\phi',cls=Function)

diffphix2 = diff(diff(phi(r),x),x)
diffphix4 = diff(diff(diffphix2,x),x)
diffphiy2 = diff(diff(phi(r),y),y)
diffphiy4 = diff(diff(diffphiy2,x),x)
diffphix2y2 = diff(diff(diffphix2,y),y)
diffphiy2x2 = diff(diff(diffphiy2,x),x)

first = diffphix2y2  +diffphiy2x2 + diffphix4 + diffphiy4
pprint(expand(first))
print '------------------'
pprint(diffphix4+diffphiy4)

wolphram = (-83*c**2*r**2*sqrt(c**2+r**r) + 6*r*(r**3*sqrt(c**2+r**2) + 20*c
#pprint(diff(diff(phi(r),x),x))
#pprint(diff(r,y))
#theta = sqrt(r**2+c**2)

#result = integrate(integrate(integrate(integrate(theta,x),x),y),y)

#pprint(result)
