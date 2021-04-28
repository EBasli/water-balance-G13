#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 08:15:49 2017

@author: theimovaara
"""

# matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import integrate

#to show you how it works 

# Definition of parameters
a = 1
b = 0.25
c = 0.1
d = 0.01


# Definition of Rate Equation
def dYdt(t, Y):
    """ Return the growth rate of fox and rabbit populations. """
    return np.array([a*Y[0] - b*Y[0]*Y[1],
                     -c*Y[1] + d*Y[0]*Y[1]])


def main():
    # Definition of output times
    tOut = np.linspace(0, 100, 200)              # time
    nOut = np.shape(tOut)[0]

    # Initial case, 10 rabbits, 5 foxes
    Y0 = np.array([10, 5])
    mt.tic()
    t_span = [tOut[0], tOut[-1]]
    YODE = sp.integrate.solve_ivp(dYdt, t_span, Y0, t_eval=tOut, 
                                  method='RK45', vectorized=True, 
                                  rtol=1e-5 )
    # infodict['message']                     # >>> 'Integration successful.'
    rODE = YODE.y[0,:]
    fODE = YODE.y[1,:]



    # Plot results with matplotlib
    plt.figure()
    plt.plot(tOut, rODE, 'r-', label='RODE')
    plt.plot(tOut, fODE, 'b-', label='FODE')


    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.ylabel('population')
    plt.title('Evolution of fox and rabbit populations')
    # f1.savefig('rabbits_and_foxes_1.png')
    plt.show()

    plt.figure()
    plt.plot(fODE, rODE, 'b-', label='ODE')


    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('Foxes')    
    plt.ylabel('Rabbits')
    plt.title('Evolution of fox and rabbit populations')
    # f2.savefig('rabbits_and_foxes_2.png')
    plt.show()


if __name__ == "__main__":
    main()
