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

r, nu = np.loadtxt("../data/8b-distrnorm.txt", comments='#', unpack=True)

#%

plt.plot(r, nu, label=r'$p_{\nu}(r)$')
plt.xlabel(r'$r = R / R_{\odot}$', fontsize=20)
plt.ylabel(r'$p_{\nu}(r)$', fontsize=20)
plt.legend(fontsize=14)
plt.title(r'${}^{8}$B Neutrino Distribution - Normalized')
plt.savefig("./nu8b-distrnorm.png", dpi=300, format='png', bbox_inches="tight")
#plt.show()

#%

print("integral =", np.trapz(nu, r))
# O calculo acima resulta em 0.9999674861119427

#%
