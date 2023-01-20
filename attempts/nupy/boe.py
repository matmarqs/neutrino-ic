#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 15:03:40 2022
@author: jaoboe
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
plt.rcParams.update({'font.size': 16})
from scipy.integrate import odeint
from scipy.interpolate import UnivariateSpline


def ODEs(y,t,fint,dfint): #y : array of functions ; t: evolution parameter

    #constants
    E = 0.1321671e6 #eV (most probable energy for a 8B solar neutrino)
    delta_m_squared = 7.42e-5 #eV^2 (difference for two neutrino approx)
    theta = 0.5836 #rad (measured mixture angle for two neutrino approx)
    G_F = 1.1663787e-23 #eV^-2



    #effective masses and matter mixing angle
    D = lambda t: 2*np.sqrt(2)*E*G_F*fint(t) #the last term is the electron density (in natural units GeV^3)
    Delta = lambda t: np.sqrt((D(t)-delta_m_squared*np.cos(2*theta))**2+(delta_m_squared*np.sin(2*theta))**2)
    theta_m_dot = lambda t: (delta_m_squared * np.sin(2*theta)*2*np.sqrt(2)*E*G_F*dfint(t)/(2*(Delta(t))))

    #functions
        #real and imaginary components of psi_1m
    psi_1m_r = y[0]
    psi_1m_i = y[1]
        #real and imaginary components of psi_2m
    psi_2m_r = y[2]
    psi_2m_i = y[3]

    #equations
        #for the components of psi_1m
    dpsi_1m_rdt = -(Delta(t)/(4*E))*psi_1m_i - theta_m_dot(t)*psi_2m_r
    dpsi_1m_idt = (Delta(t)/(4*E))*psi_1m_r - theta_m_dot(t)*psi_2m_i

        #for the components of psi_2m
    dpsi_2m_rdt = theta_m_dot(t)*psi_1m_r + (Delta(t)/(4*E))*psi_2m_i
    dpsi_2m_idt = theta_m_dot(t)*psi_1m_i - (Delta(t)/(4*E))*psi_2m_r

    return [dpsi_1m_rdt,dpsi_1m_idt,dpsi_2m_rdt,dpsi_2m_idt]

#solar electron density
data = pd.read_csv("/home/jaoboe/Desktop/programas-ic/programas-python/dados-densidade-solar/dados-densidade-solar.csv",
                   names = ["t","n_e"]) #data from John Bahcall

t = np.float64(data["t"]) #actually radius in solar radius units
ne =10**(27)* np.float64(data["n_e"]) #number of electron * eV^3

n_e = UnivariateSpline(t, ne,k = 3,s=0)
dn_e = n_e.derivative(n=1)

#set initial conditions
y0 = [1/np.sqrt(2),0,1/np.sqrt(2),0]

#set time range (still need to actually evaluate the time scale for this problem)
t = np.linspace(0,1,28)

#integrating the system of diff equations
y = odeint(ODEs,y0,t,args = (n_e,dn_e))

psi_1m_r = y[:,0]
psi_1m_i = y[:,1]

psi_2m_r = y[:,2]
psi_2m_i = y[:,3]

print(psi_1m_r**2+psi_1m_i**2)
fig1,ax1 = plt.subplots(figsize = (10,6))
ax1.plot(t,psi_1m_r**2+psi_1m_i**2,label = r"$|\psi_{1m}|^2$")
ax1.plot(t,psi_2m_r**2+psi_2m_i**2,label = r"$|\psi_{2m}|^2$")
ax1.set_xlabel(r"$R/R_\odot$")
ax1.set_ylabel(r"$|\psi|^2$")
ax1.legend()

fig1.savefig("adiabatic-condition-two-generation.png")
plt.show()
