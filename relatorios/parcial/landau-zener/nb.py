#%

import numpy as np
import matplotlib.pyplot as plt
from sympy import *
from matplotlib import rcParams
rcParams['text.usetex'] = True
plt.style.use('bmh')
init_printing()

#%

### Funcoes arbitrarias
x = symbols('x')
f, g, m = symbols('f g m', cls = Function)
f = 30 + 2*(x-3)**2 + 10*x
g = 20*exp(-((x-3)/2)**2) + 10*x
fx = lambdify(x, f, modules = ['numpy'])
gx = lambdify(x, g, modules = ['numpy'])
xvar = np.linspace(-3, 9, 1000)
plt.plot(xvar, fx(xvar))
plt.plot(xvar, gx(xvar))
plt.show()

#%

### Retirando a media das funcoes
x = symbols('x')
f, g, m, d = symbols('f g m d', cls = Function)
f = 30 + 2*(x-3)**2 + 10*x
g = 20*exp(-((x-3)/2)**2) + 10*x
m = (f + g)/2
fx = lambdify(x, f, modules = ['numpy'])
gx = lambdify(x, g, modules = ['numpy'])
mx = lambdify(x, m, modules = ['numpy'])
xvar = np.linspace(-9 - 3, 15 - 3, 1000)
plt.plot(xvar, fx(xvar + 3) - mx(xvar + 3))
plt.plot(xvar, gx(xvar + 3) - mx(xvar + 3))
plt.show()

#%

### Diferenca das funcoes
# f* = (f(x) - g(x))/2
# g* = (g(x) - f(x))/2
# f* - g* = f(x) - g(x)
d = f - g
dx = lambdify(x, d, modules = ['numpy'])
xvar2 = np.linspace(-3, 8, 1000)
plt.plot(xvar2, dx(xvar2))
plt.show()

#%

### Energias do problema do Landau-Zener
x = symbols('x')
a, D = 1, 1
E1, E2 = symbols('E1 E2')
E1 = (1/2)*sqrt( (a*x)**2 + D**2 )
E2 = -(1/2)*sqrt( (a*x)**2 + D**2 )
E1x = lambdify(x, E1, modules = ['numpy'])
E2x = lambdify(x, E2, modules = ['numpy'])
xvarE = np.linspace(-8, 8, 1000)
plt.plot(xvarE,E1x(xvarE),label=r'$E_+$')
plt.plot(xvarE,E2x(xvarE),label=r'$E_-$')
plt.plot(xvarE,  (1/2)*a*xvarE, color='#444444', linestyle='--', dashes=(5,6), linewidth=1.5)
plt.plot(xvarE, -(1/2)*a*xvarE, color='#444444', linestyle='--', dashes=(5,6), linewidth=1.5)
plt.ylabel(r'$E$', fontsize=20)
plt.xlabel(r'$t$', fontsize=20)
plt.legend(fontsize=14)
plt.savefig('figures/gap.png', dpi=300, format='png', bbox_inches="tight")
plt.show()

#%

dE = E1 - E2
dEx = lambdify(x, dE, modules = ['numpy'])
xvar2 = np.linspace(-3, 8, 1000)
plt.plot(xvarE, dEx(xvarE))
plt.show()

#%

xlin = np.linspace(-5, 5, 500)
plt.plot(xlin,  xlin, label=r'$E_+ = +\frac{\alpha t}{2}$')
plt.plot(xlin, -xlin, label=r'$E_- = -\frac{\alpha t}{2}$')
plt.ylabel(r'$E$', fontsize=20)
plt.xlabel(r'$t$', fontsize=20)
plt.legend(fontsize=14)
plt.savefig('figures/nogap.png', dpi=300, format='png', bbox_inches="tight")
plt.show()

#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%


#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%





#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%


#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%


#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%


#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%


#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%


#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%



#%


