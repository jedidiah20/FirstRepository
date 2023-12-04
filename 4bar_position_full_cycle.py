# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 20:53:59 2020

@author: Adam Wickenheiser
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root_scalar

# givens
R1 = 8                               # ground link length [cm]
R2 = 3.75                            # crank length [cm]
R3 = 8.5                             # coupler length [cm]
R4 = 5.25                            # rocker length [cm]
w = np.pi                            # rotation rate [rad/s]
n = 1                                # number of cycles to plot
link_config = 'open'                 # 'open' or 'crossed'
z1 = R1                              # ground link

# set up arrays for simulation
t = np.linspace(0,n*2*np.pi/w,101)   # time array [s]
tht2 = np.zeros_like(t)              # crank angle array [rad]
tht3 = np.zeros_like(t)              # coupler angle array [rad]
tht4 = np.zeros_like(t)              # rocker angle array [rad]


# function whose root we want to find within the bracket
def calc_tht4(tht4_guess):
    z4 = R4*np.exp(1j*tht4_guess)
    z5 = z1 + z4
    z3 = z5 - z2
    r = np.abs(z3)
    return r - R3


# loop through all time steps and calculate angles
for i in range(t.size):
    tht2[i] = w*t[i]
    z2 = R2*np.exp(1j*tht2[i])

    # choose which solution of tht4 to look for
    tht_div = np.angle(z2-z1)
    if link_config == 'open':
        if tht2[i] >= 0:
            bracket = [tht_div-np.pi,tht_div]      # angle range to search within
            x0 = tht_div-np.pi/2                   # initial guess for tht4
        else:
            bracket = [tht_div,tht_div+np.pi]      # angle range to search within
            x0 = tht_div+np.pi/2                   # initial guess for tht4
    else:  # 'crossed'
        if tht2[i] >= 0:
            bracket = [tht_div,tht_div+np.pi]      # angle range to search within
            x0 = tht_div+np.pi/2                   # initial guess for tht4
        else:
            bracket = [tht_div-np.pi,tht_div]      # angle range to search within
            x0 = tht_div-np.pi/2                   # initial guess for tht4

    sol = root_scalar(calc_tht4,x0=x0,bracket=bracket)
    tht4[i] = sol.root
    z4 = R4*np.exp(1j*tht4[i])
    z5 = z1 + z4
    z3 = z5 - z2
    tht3[i] = np.angle(z3)

# positions of A and B
z_A = R2*np.exp(1j*tht2)
z_B = z_A + R3*np.exp(1j*tht3)

plt.figure()
plt.plot(np.real(z_A),np.imag(z_A),label='Point A')
plt.plot(np.real(z_B),np.imag(z_B),label='Point B')
plt.axis('equal')
plt.legend()

