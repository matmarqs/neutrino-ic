#%

#%matplotlib
import numpy as np
import matplotlib.pyplot as plt
#from sympy import *
#import sympy as sp
from matplotlib import rcParams
rcParams['text.usetex'] = True
plt.style.use('bmh')
#init_printing()

#%

E, nu, nu_sup, nu_inf = np.loadtxt("../data/8b-energy.txt", comments='#', unpack=True)

#%

plt.plot(E, nu, label=r'$\nu(E)$')
plt.xlabel(r'$E$', fontsize=20)
plt.ylabel(r'$\nu(E)$', fontsize=20)
plt.legend(fontsize=14)
plt.title(r'${}^{8}$B Neutrino Energy Spectrum')
plt.savefig("./nu8b-energy.png", dpi=300, format='png', bbox_inches="tight")
#plt.show()

#%

print("integral =", np.trapz(nu, E))
# O calculo acima resulta em 0.999986544

#%
