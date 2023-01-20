rang = 60000
th_graus = 33.44
dm2 = 7.42e-5       # NuFit
factor = 1.2669326791371    # 1.27 da conversão de unidades

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
plt.style.use('bmh')

#rcParams['text.usetex'] = True
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def graustorad(a):
    return (a * np.pi) / 180

# vacuum constants
th = graustorad(th_graus)  # theta12, do NuFit

sin2th = (np.sin(2*th))**2
Delta = dm2
#Delta = np.sqrt( (D - dm2 * np.cos(2*th))**2 + (dm2 * np.sin(2*th))**2 )

def P_e(LE):    # LE significa L/E
    return 1 - sin2th  * np.sin(factor * dm2 * LE)**2

LE_array = np.linspace(0, rang, 500)    # L/E em km/GeV

plt.plot(LE_array , P_e(LE_array), label=r'$P(\nu_e \to \nu_e)$')
plt.plot(LE_array , 1-P_e(LE_array), label=r'$P(\nu_e \to \nu_\mu)$')
plt.xlabel(r'$L/E$ (km/GeV)', fontsize=20)
plt.ylabel(r'Probabilidade', fontsize=20)
plt.ylim(0, 1.05)
plt.legend(fontsize=16)
plt.title(r"Oscilação no vácuo para 2 gerações de neutrinos")
plt.savefig("./2nu-vacuo.png", dpi=300, format='png', bbox_inches="tight")
