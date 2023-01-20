#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys, os
import numpy as np
import pandas as pd
from scipy.integrate import quad, dblquad
from matplotlib import pyplot as plt
from matplotlib import rc
plt.style.use('bmh')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


# In[2]:


energ = np.loadtxt('events.txt', comments='#')
kamland = pd.DataFrame(data=sorted(energ), columns=['nu energy (MeV)'])
dbin = ( energ.max() - energ.min() ) / 13


#plt.hist(energ, bins=13)
#plt.show()


# data taking from March 9, 2002 to January 11, 2004
# 2002 é bissexto
time_0 = 257.3   # days
time_1 = 257.9   # days
life_time = time_0 + time_1
n_p = 4.61e+31      # number of target photons in the fiducial volume
eff = 89.8 / 100
sys_desvpad = 6.5 / 100
N_T = n_p * eff * life_time     # normalization constant
resol = 6.75 / 100

m_p = 938.27208816  # MeV / c^2
m_n = 939.56542052  # MeV / c^2
m_e = 0.51099895000 # MeV / c^2
dM = m_n - m_p - m_e


def sig(e):
    return resol / np.sqrt(e - dM)

def R(e, ei):   # inputs e, ei in MeV
    sg = sig(e)
    return (1 / (sg * np.sqrt(2*np.pi))) * np.exp(-(e-ei)**2/(2*sg**2))


#Ei_min = energ.min() + dM
#Ei_max = energ.max() + dM
#u = np.linspace(Ei_min, Ei_max, 1000)
#v = np.linspace(Ei_min, Ei_max, 1000)
#U, V = np.meshgrid(u, v)
#W = R(U, V)
#fig = plt.figure()
#ax = fig.add_subplot(projection='3d')
#ax.plot_surface(U, V, W)
#plt.show()


ordering = 'normal'
# parameters from NuFIT
th12 = rad(33.44)
th23 = rad(49.2)
th13 = rad(8.57)
d_cp = rad(194.0)
dm2_21 = 7.42e-5    # eV^{-2}
dm2_atm = 2.515e-3  # eV^{-2}
if ordering == 'inverted':
    dm2_32 = dm2_atm
else:
    dm2_32 = - dm2_atm - dm2_21
dm2_13 = - dm2_21 - dm2_32

def sci(num):
    return "{:e}".format(num)

def rad(graus):
    return np.pi * graus / 180.0


def P_e0(E, L):
    return 1.0


def P_e(E, L):
    a = 1266.9326791370845 * dm2_21 * L / E
    return 1 - (np.sin(2*th12) * np.sin(a))**2


#def P_e(E, L):
#    d13 = 1266.9326791370845 * dm2_13 * L / E   # 1266.93 is the conversion factor for dm^2 L / 4E
#    d21 = 1266.9326791370845 * dm2_21 * L / E   # E (MeV), L (km), dm2 (eV^2)
#    d32 = 1266.9326791370845 * dm2_32 * L / E
#    p = np.cos(th12)**2 * np.sin(d13)**2 + np.sin(th12)**2 * np.sin(d32)**2
#    p *= - np.sin(2*th13)**2
#    p_o = np.sin(2*th12)**2 * np.cos(th13)**4 * np.sin(d21)**2
#    p += 1 - p_o
#    return p


def cross_sec(E):   # E in MeV
    E_e = E - dM
    p_e = np.sqrt(E_e**2 - m_e**2)
    return 1e-1 * p_e * E_e * E**(-0.07056 + 0.02018*np.log(E)-0.001953*(np.log(E))**3)  # return in 1e+42 cm^2

def cross_sec0(E):   # E in MeV
    E_e = E - dM
    p_e = np.sqrt(E_e**2 - m_e**2)
    return 0.0952 * E_e * p_e  # return in 1e+42 cm^2 (have to multiply it by 1e-42)


#x = np.linspace(dM + m_e, 10, 200)
#plt.plot(x, cross_sec0(x), label='naive')
#plt.plot(x, cross_sec(x), label='improved')
#plt.legend()
#plt.show()


isotop = pd.DataFrame({'235U':[3.217, -3.111, 1.395, -3.690e-1, 4.445e-2, -2.053e-3],
        '238U':[4.833e-1, 1.927e-1, -1.283e-1, -6.762e-3, 2.233e-3, -1.536e-4],
        '239Pu':[6.413, -7.432, 3.535, -8.820e-1, 1.025e-1, -4.550e-3],
        '241Pu':[3.251, -3.204, 1.428, -3.675e-1, 4.254e-2, -1.896e-3]})
isotop.loc['Q']    = [202.36, 205.99, 211.12, 214.26]  # 0.26, 0.52, 0.34, 0.33 (incertezas)
isotop.loc['PWR']  = [ 0.560,  0.080,  0.300, 0.060 ]
isotop.loc['MOX']  = [ 0.000,  0.081,  0.708, 0.212 ]
isotop.loc['PHWR'] = [ 0.543,  0.411,  0.022, 0.024 ]
lambd = {}
lambd['235U']  = lambda z: np.exp(sum([isotop['235U'][j] * z**j for j in range(6)]))
lambd['238U']  = lambda z: np.exp(sum([isotop['238U'][j] * z**j for j in range(6)]))
lambd['239Pu'] = lambda z: np.exp(sum([isotop['239Pu'][j] * z**j for j in range(6)]))
lambd['241Pu'] = lambda z: np.exp(sum([isotop['241Pu'][j] * z**j for j in range(6)]))


#X = np.linspace(2, 8, 200)
#for isot in isotop:
#    plt.plot(X, lambd[isot](X), label=isot)
#    plt.legend()
#plt.yscale('log')
#plt.show()


class Reactor:
    def __init__(self, nm, typ, pth, d, LF_2002, LF_2003):
        self.name = nm
        self.type = typ
        self.P_th = pth
        self.dist = d
        self.LF_0 = LF_2002
        self.LF_1 = LF_2003
    def df(self):
        c = ['type', 'P_th (MW)', 'dist (km)', 'LF_0 (%)', 'LF_1 (%)']
        dataf = pd.DataFrame(columns=c)
        dataf.loc[self.name] = [self.type, self.P_th, self.dist,
                           self.LF_0, self.LF_1]
        return dataf
    def LF(self):
        lf = (self.LF_0 * time_0 + self.LF_1 * time_1) / (time_0 + time_1)
        return lf / 100  # because LF was in %
    def p(self, i):
        pwr = {'235U' : 0.560,
               '238U' : 0.080,
               '239Pu': 0.300,
               '241Pu': 0.060, }
        mox = {'235U' : 0.000,
               '238U' : 0.081,
               '239Pu': 0.708,
               '241Pu': 0.212, }
        phwr = {'235U' : 0.543,
                '238U' : 0.411,
                '239Pu': 0.022,
                '241Pu': 0.024, }
        if self.type == 'PHWR':
            return phwr[i]
        elif self.type == 'MOX':
            return mox[i]
        else:
            return pwr[i]


def create_reactor(file_name):
    f = open(file_name, 'r')
    name = f.readline().split()[1]
    typ  = f.readline().split()[1]
    pth  = float(f.readline().split()[1])
    d    = float(f.readline().split()[1])
    f.readline()
    for line in f:
        year = int(line.split()[0])
        if year == 2002:
            break
    LF_2002 = float(line.split()[7])
    LF_2003 = float(f.readline().split()[7])
    f.close()
    return Reactor(name, typ, pth, d, LF_2002, LF_2003)


reactors = []
for file in os.listdir('../3-data/reactors/'):
    if file.endswith('.react'):
        r = create_reactor('../3-data/reactors/' + file)
        reactors.append(r)
reac = pd.concat([r.df() for r in reactors])
print('reactors dataframe shape =', reac.shape)


def phi(E):
    soma = 0
    for r in reactors:
        a = P_e(E, r.dist) * r.P_th * r.LF() / (4 * np.pi * r.dist**2)
        b = 0
        for i in isotop:
            b += r.p(i) * lambd[i](E) / isotop.loc['Q'][i]
        soma += a * b
    return soma

def phi0(E):
    soma = 0
    for r in reactors:
        a = P_e0(E, r.dist) * r.P_th * r.LF() / (4 * np.pi * r.dist**2)
        b = 0
        for i in isotop:
            b += r.p(i) * lambd[i](E) / isotop.loc['Q'][i]
        soma += a * b
    return soma


days = 24 * 60 * 60 # in seconds
elec_charge = 1.60217663e-19  # electron charge in coulombs
M = 1e-42 * N_T * 1e-10 * days / elec_charge  # constant to make units right

def func(E_nu, Ei):  # f(y, x)
     return M * phi(E_nu) * cross_sec(E_nu) * R(E_nu, Ei)
def func0(E_nu, Ei):  # f(y, x)
     return M * phi0(E_nu) * cross_sec(E_nu) * R(E_nu, Ei)


Ei_min = energ.min() + dM
Ei_max = energ.max() + dM
print('Ei_min =', Ei_min)
print('Ei_max =', Ei_max)

# integrating on [-6 sigma, 6 sigma] because that is sufficient and faster
N0 = dblquad(func0, Ei_min, Ei_max, lambda x: x-6*sig(x), lambda x: x+6*sig(x))
print('N0 =', N0)

N = dblquad(func, Ei_min, Ei_max, lambda x: x-6*sig(x), lambda x: x+6*sig(x))
print('N =', N)


print(N[0], '+-', 0.08 * N[0])
print(N0[0], '+-', 0.08 * N0[0])
print(N[0] / N0[0], '+-', 0.08 * N[0]/N0[0])
