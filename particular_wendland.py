from __future__ import division
from sympy import *
from sympy.plotting import plot

from sympy.functions import exp
from sympy.core.containers import Tuple

x,y,rs = symbols('x y r',real=True)
c = symbols('c')
r = sqrt(x**2+y**2)
pprint(diff(r,x))
pprint(diff(r,y))
theta = sqrt(r**2+c**2)

#result = integrate(integrate(integrate(integrate(theta,x),x),y),y)

#pprint(result)
