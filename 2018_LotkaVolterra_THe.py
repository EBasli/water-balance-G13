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
import MyTicToc as mt

def tic():
    # Homemade version of matlab tic and toc functions
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()


def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print("Elapsed time is " + str(time.time() - startTime_for_tictoc)
              + " seconds.")
    else:
        print("Toc: start time not set")

#Characteristics Cell VP-06 Wieringermeer landfill:
#base area [m2]
A_base = 28355
#top area [m2]
A_top = 9100
#slope witdh [m]
w = 38
#waste body height [m]
h_body = 12
#cover layer height [m]
h_cover = 1.5
#waste (net weight [kg])
W = 281083000


#Data import
import pandas as pd

meteo = pd.read_excel('WieringermeerData_Meteo.xlsx')
leachate = pd.read_excel('WieringermeerData_LeachateProduction.xlsx')

# Definition of data
rain = meteo.rain_station       #m/day
evap = meteo.pEV                #m/day
temp = meteo.temp               #celsius
leach_cum = leachate.loc[:,0]   #m3/day   CUMULATIVE
leach = []                      #m3/day
leach.append(leach_cum[0])
for i in range(len(leach_cum)-1):
    l = leach_cum[i+1] - leach_cum[i]
    leach.append(l)



# Definition of Rate Equation
def dSdt(t, S):
    return np.array([S[0] - S[2] - S[1],
                     (1-S[4])*S[2] - S[3],
                     S[4]*S[2] + S[3] - S[5]])
#S = [J_rf, E, L_cl, L_wb, beta, Q_dr]




def water_balance():
    #Definition of time
    tOut = np.linspace(0, 100, 200)
    nOut = np.shape(tOut)[0]
    
    #Initial conditions
    S_cl = 10
    S_wb = 10
    mt.tic()
    t_span = [tOut[0], tOut[-1]]
    
    #Definitions
    a = 10  # saturated hydraulic conductivity (m/day)
    S_cl_min = 10  # min achievable storage in cover layer
    S_cl_max = 10  # max achievable storage in cover layer
    S_wb_min = 10  # min achievable storage in waste layer
    S_wb_max = 10  # max achievable storage in waste layer
    b_cl = 10  # emprical parameter (dimensionless)
    b_wb = 10  # emprical parameter (dimensionless)
    pEv = evap # potential evaporation (m/day) from excel file
    C_f = 10     # crop factor
    S_EV_min = 10
    S_EV_max = 10
    beta_0 = 10    # beta zero
    
    #Equations
    J_rf = rain
    if S_cl < S_EV_min:
        f_red = 0
    elif S_cl > S_EV_max:
        f_red = 1
    else:
        (S_cl - S_EV_min) / (S_EV_max - S_EV_min)
    E = pEv * C_f * f_red
    L_cl = a * ((S_cl - S_cl_min) / (S_cl_max - S_cl_min))**b_cl
    L_wb = a * ((S_wb - S_wb_min) / (S_wb_max - S_wb_min))**b_wb
    beta = beta_0 * ((S_cl - S_cl_min) / (S_cl_max - S_cl_min))
    Q_dr = beta * L_cl + L_wb
    S = [J_rf, E, L_cl, L_wb, beta, Q_dr]
    
    
    

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

    '''EULER'''
    # Initialize output vector
    YEuler = np.zeros([nOut, 2], dtype=float)

    dtMax = 0.1
    # dtMin = 1e-11
    t = tOut[0]
    iiOut = 0

    # Initialize problem
    mt.tic()
    Y = Y0
    # Write initial values to output vector
    YEuler[iiOut, [0, 1]] = Y
    while (t < tOut[nOut-1]):
        # check time steps
        Rates = dYdt(t, Y)
        dtRate = -0.7 * Y/(Rates + 1e-18)
        dtOut = tOut[iiOut+1]-t
        dtRate = (dtRate <= 0)*dtMax + (dtRate > 0)*dtRate
        dt = min(min(dtRate), dtOut, dtMax)

        Y = Y + Rates * dt
        t = t + dt

        # print ("Time to Output is " + str(np.abs(tOut[iiOut+1]-t)) +
        # " days.")

        if (np.abs(tOut[iiOut+1]-t) < 1e-5):
            YEuler[iiOut+1, [0, 1]] = Y
            iiOut += 1

    rEuler, fEuler = YEuler.T
    mt.toc()

    '''EULER Predictor Corrector'''
    # Initialize output vector
    YPC = np.zeros([nOut, 2], dtype=float)

    dtMax = 0.01
    # dtMin = 1e-11
    t = tOut[0]
    iiOut = 0

    # Initialize problem
    mt.tic()
    Y = Y0
    # Write initial values to output vector
    YPC[iiOut, [0, 1]] = Y
    while (t < tOut[nOut-1]):
        # check time steps
        Rates = dYdt(t, Y)
        dtRate = -0.7 * Y/(Rates + 1e-18)
        # We only need to take the positive rates in to account
        dtRate = (dtRate < 0)*dtMax + (dtRate >= 0)*dtRate
        dtOut = tOut[iiOut+1]-t
        dt = min(min(dtRate), dtOut, dtMax)

        # Iteration with predictor corrector
        converged = False
        convCrit = 1e-6
        Yn = Y + Rates * dt
        while not(converged):
            Rn = dYdt(t, Yn)
            Ynn = Y + (Rates + Rn)/2 * dt
            if max(np.abs(Yn-Ynn)) < convCrit:
                converged = True
            else:
                Yn = Ynn

        Y = Yn
        t = t + dt
    # print ("Time to Output is " + str(np.abs(tOut[iiOut+1]-t)) + " days.")

        if (np.abs(tOut[iiOut+1]-t) < 1e-5):
            YPC[iiOut+1, [0, 1]] = Y
            iiOut += 1

    rPC, fPC = YPC.T
    mt.toc()

    '''RungeKutta'''
    # Initialize output vector
    YRK = np.zeros([nOut, 2], dtype=float)
    dtMax = 0.01
    # dtMin = 1e-11
    t = tOut[0]
    iiOut = 0
    # Initialize problem
    mt.tic()
    Y = Y0
    # Write initial values to output vector
    YRK[iiOut, [0, 1]] = Y
    while (t < tOut[nOut-1]):
        # check time steps
        R1 = dYdt(t, Y)

        dtRate = -0.7 * Y/(R1 + 1e-18)
        # We only need to take the positive rates in to account
        dtRate = (dtRate < 0)*dtMax + (dtRate >= 0)*dtRate
        dtOut = tOut[iiOut+1]-t
        dt = min(min(dtRate), dtOut, dtMax)

        k1 = dt*R1
        k2 = dt*dYdt(t+dt/2, Y+k1/2)
        k3 = dt*dYdt(t+dt/2, Y+k2/2)
        k4 = dt*dYdt(t+dt, Y+k3)

        Y = Y+(k1+2*k2+2*k3+k4) / 6
        t = t + dt

        # print ("Time to Output is " + str(np.abs(tOut[iiOut+1]-t))
        # + " days.")
        if (np.abs(tOut[iiOut+1]-t) < 1e-5):
            YRK[iiOut+1, [0, 1]] = Y
            iiOut += 1

    rRK, fRK = YRK.T
    mt.toc()

    # Plot results with matplotlib
    plt.figure()
    plt.plot(tOut, rODE, 'r-', label='RODE')
    plt.plot(tOut, fODE, 'b-', label='FODE')
    plt.plot(tOut, rEuler, 'g+', label='REuler')
    plt.plot(tOut, fEuler, 'm+', label='FEuler')
    plt.plot(tOut, rPC, 'rx', label='RPC')
    plt.plot(tOut, fPC, 'bx', label='FPC')
    plt.plot(tOut, rRK, 'g.', label='RRK')
    plt.plot(tOut, fRK, 'm.', label='FRK')

    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.ylabel('population')
    plt.title('Evolution of fox and rabbit populations')
    # f1.savefig('rabbits_and_foxes_1.png')
    plt.show()

    plt.figure()
    plt.plot(fODE, rODE, 'b-', label='ODE')
    plt.plot(fEuler, rEuler, 'b+', label='Euler')
    plt.plot(fPC, rPC, 'r-', label='Predictor Corrector')
    plt.plot(fRK, rRK, 'g-', label='Runge Kutta')

    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('Foxes')    
    plt.ylabel('Rabbits')
    plt.title('Evolution of fox and rabbit populations')
    # f2.savefig('rabbits_and_foxes_2.png')
    plt.show()


if __name__ == "__main__":
    main()
