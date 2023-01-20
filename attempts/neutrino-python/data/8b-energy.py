import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

plt.style.use('bmh')

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


E, nu, nu_sup, nu_inf = np.loadtxt("../data/8b-energy.txt", comments='#', unpack=True)


plt.plot(E, nu, label=r'$p_{E}(E)$')
plt.xlabel(r'$E$', fontsize=20)
plt.ylabel(r'$p_{E}(E)$', fontsize=20)
plt.legend(fontsize=16)
plt.title(r'Densidade de probabilidade de produção com energia $E$ para neutrinos $^8$B')
plt.savefig("./nu8b-energy.png", dpi=300, format='png', bbox_inches="tight")

print("integral =", np.trapz(nu, E))
# O calculo acima resulta em 0.999986544
