#%

#%matplotlib
import numpy as np
#import matplotlib.pyplot as plt
#from sympy import *
#import sympy as sp
#from matplotlib import rcParams
#rcParams['text.usetex'] = True
#plt.style.use('bmh')
#init_printing()

#from scipy import interpolate
#from scipy import integrate

#%

r, nu = np.loadtxt("../data/8b-distr.txt", comments='#', unpack=True)

#f = interpolate.interp1d(r, nu, kind = 'cubic')

#%

#integral, abserr = integrate.quad(f, 0.0, r[len(r)-1])
#print("integral = {:.8e}, abserr = {:.8e}".format(integral, abserr))
# result:
# integral = 4.09847562e-04, abserr = 9.72090834e-09

normalization = 4.09847562e-04

Ndata = len(r)
for i in range(Ndata):
    print(" {:.5f}   {:.10e}".format(r[i], nu[i] / normalization))

#print("f({:.7f}) = {:.7e}".format(-0.00082, f(-0.00082)))
#print("f({:.7f}) = {:.7e}".format(-0.00041, f(-0.00041)))
#print("f({:.7f}) = {:.7e}".format(0, f(0)))
#print("f({:.7f}) = {:.7e}".format(0.00041, f(0.00041)))
#print("f({:.7f}) = {:.7e}".format(0.00082, f(0.00082)))

#plt.plot(r, f(r), label=r'$p_{\nu}(r)$')
#plt.xlabel(r'$r = R / R_{\odot}$', fontsize=20)
#plt.ylabel(r'$p_{\nu}(r)$', fontsize=20)
#plt.legend(fontsize=14)
#plt.title(r'${}^{8}$B Neutrino Distribution')
#plt.savefig("./nu-interp.png", dpi=300, format='png', bbox_inches="tight")
#plt.show()

#%

#print("integral =", np.trapz(nu, r))
# O calculo acima resulta em 0.00040983303311751916

#%
