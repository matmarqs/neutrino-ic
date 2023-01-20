import numpy as np
import matplotlib.pyplot as plt

r, elecdens = np.loadtxt("./elecdens.txt", comments='#', unpack=True)

media = sum(elecdens)/len(elecdens)

for _ in r:
    print(("{:.7f}".format(_))+("   {:.6f}".format(media)))

plt.plot(r, elecdens, label='real data')
plt.plot(r, np.full(len(r), media), label='media')
plt.xlabel(r'$R/R_{\odot}$', fontsize=20)
plt.ylabel(r'$N_e/N_A$', fontsize=20)
plt.legend(fontsize=14)
plt.title(r'Electron density')
plt.savefig("elecdens-media.png", dpi=300, format='png', bbox_inches="tight")
