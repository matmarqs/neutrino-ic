import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

plt.style.use('bmh')

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


r, nu = np.loadtxt("../data/8b-distrnorm.txt", comments='#', unpack=True)


plt.plot(r, nu, label=r'$p_{r}(r)$', color='brown')
plt.xlabel(r'$r = R / R_{\odot}$', fontsize=20)
plt.ylabel(r'$p_{r}(r)$', fontsize=20)
plt.legend(fontsize=16)
plt.title(r'Densidade de probabilidade de produção em $r = R/R_{\odot}$ para neutrinos $^8$B')
plt.savefig("./nu8b-raio.png", dpi=300, format='png', bbox_inches="tight")
#plt.show()


print("integral =", np.trapz(nu, r))
# O calculo acima resulta em 0.9999674861119427
