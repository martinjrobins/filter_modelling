from sympy import *

#2D1P greens function

x,y,kx,ky,mu = symbols('x y kx ky mu')
r = sqrt(x**2+y**2)

A = 0.5*log(2*(cosh(ky)-cos(kx)))
Ax = simplify(diff(A,kx))
Ay = simplify(diff(A,ky))
Axx = simplify(diff(Ax,kx))
Axy = simplify(diff(Ax,ky))

print "A = "
pprint(A)
print "Ax = "
pprint(Ax)
print "Ay = "
pprint(Ay)
print "Axx = "
pprint(Axx)
print "Axy = "
pprint(Axy)

Sxx = simplify(A + ky*Ay - 1)
SxxST = simplify(log(r) - x*x/r**2)
Sxy = simplify(-ky*Ax)
Syy = simplify(A - ky*Ay)
SyyST = simplify(log(r) - x*x/r**2)
Tyxy = 2*mu*simplify(Ax - (ky*Axy+Ax))
Tyyy = 2*mu*simplify(ky*Axx + Ay)

print "Sxx = "
pprint(Sxx)
print "Sxy = "
pprint(Sxy)
print "Syy = "
pprint(Syy)
print "Tyxy = "
pprint(Tyxy)
print "Tyyy = "
pprint(Tyyy)

#Analytical integrals for Stokelet singularity
dx,dy,h,t = symbols('d_x d_y h t')

x = dx*h*t
y = dy*h*t
r = h*t
Sxx = simplify(       - x*x/r**2)
Sxy = simplify(       - x*y/r**2)
Syy = simplify(       - y*y/r**2)
Ixx = simplify(2*integrate(log(r),(t,0,0.5)) + integrate(Sxx,(t,-0.5,0.5)))
Ixy = simplify(integrate(Sxy,(t,-0.5,0.5)))
Iyy = simplify(2*integrate(log(r),(t,0,0.5)) + integrate(Syy,(t,-0.5,0.5)))
Tyxy = 4.0*mu*simplify(y*x*y/r**4)
Tyyy = 4.0*mu*simplify(y*y*y/r**4)
Iyxy = simplify(integrate(Tyxy,(t,-0.5,0.5)))
Iyyy = simplify(integrate(Tyyy,(t,-0.5,0.5)))


print "Sxx = "
pprint(Sxx)
print "Ixx = "
pprint(Ixx)
print "Ixy = "
pprint(Ixy)
print "Iyy = "
pprint(Iyy)
print "Iyxy = "
pprint(Iyxy)
print "Iyyy = "
pprint(Iyyy)

print "Tyxy = "
pprint(Tyxy)
print "Tyyy = "
pprint(Tyyy)


#Analytical integrals for 2D laplace singularity
dx,dy,h,t,k = symbols('d_x d_y h t k')

x = dx*h*t
y = dy*h*t
A = 0.5*log(2*(cosh(k*y)-cos(k*x)))
Ax = -0.5*k*sin(k*x)/(cos(k*x)-cosh(k*y))
Ay = -0.5*k*sinh(k*y)/(cos(k*x)-cosh(k*y))
IAx = simplify(integrate(Ax,(t,-0.5,0.5)))
IAy = simplify(integrate(Ay,(t,-0.5,0.5)))

print "IAx = "
pprint(IAx)
print "IAy = "
pprint(IAy)


