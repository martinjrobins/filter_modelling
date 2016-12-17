from __future__ import division
from sympy import *
from sympy.plotting import plot

from sympy.functions import exp
from sympy.core.containers import Tuple


xi = symbols(r'x_i',real=True)
yi = symbols(r'y_i',real=True)

xj = symbols(r'x_j',real=True)
yj = symbols(r'y_j',real=True)
ri = Matrix([xi,yi])
rj = Matrix([xj,yj])
dx = ri-rj



hi = symbols(r'h_i',positive=true)
hj = symbols(r'h_j',positive=true)
h = 0.5*(hi+hj)

d = symbols(r'd',positive=true)
k = symbols(r'k',positive=true)

#kij = exp(-dx.dot(dx)/(h**2))
kij = Piecewise(((2-dx.norm()/h)**4*(1 + 2*dx.norm()/h),dx.norm()<2*h),(0,True))

print '---------------------------------------'
print '           k_ij'
print '---------------------------------------'
pprint(kij)
pprint(simplify(kij))

laplace = diff(diff(kij,xi),xi) + diff(diff(kij,yi),yi)

print '-----------------------------------------------------'
print '           [dfdx,dfdy]'
print '----------------------------------------------------'

gradient = Matrix([diff(kij,xi),diff(kij,yi)])
pprint(simplify(gradient))

print '-----------------------------------------------------'
print '           [df2dx2,df2dxy]'
print '           [df2dyx,df2dy2]'
print '----------------------------------------------------'
laplace = Matrix([[diff(diff(kij,xi),xi),diff(diff(kij,yi),xi)],
                  [diff(diff(kij,xi),yi),diff(diff(kij,yi),yi)]])

print '-----------------------------------------------------'
print '           df2dx2'
print '----------------------------------------------------'
pprint(simplify(laplace[0,0]))
print 'limit is '
pprint(4*80/(hi + hj)**2)
small = 0.0000000001
#lim= laplace[0,0].subs(xi,0.0).subs(yi,0.0).subs(yj,small).subs(xj,small)
#pprint(simplify(expand(lim)))

print '-----------------------------------------------------'
print '           df2dy2'
print '----------------------------------------------------'
pprint(simplify(laplace[1,1]))
print 'limit is '
pprint(4*80/(hi + hj)**2)
#lim= laplace[1,1].subs(xi,0.0).subs(yi,0.0).subs(yj,small).subs(xj,small)
#pprint(simplify(expand(lim)))

print '-----------------------------------------------------'
print '           dfdxy'
print '----------------------------------------------------'
pprint(simplify(laplace[0,1]))
print 'limit is 0'

