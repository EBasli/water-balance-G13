# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 14:15:26 2021

@author: theimovaara
"""



# matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import integrate

# Definition of parameters
a = 1
b = 0.25
c = 0.1
d = 0.01


# Definition of Rate Equation
def dYdt(t, Y):
    """ Return the growth rate of fox and rabbit populations. """
    return np.array([ a * Y[0] - b * Y[0] * Y[1],
                     -c * Y[1] + d * Y[0] * Y[1]])


def main():
    # Definition of output times
    tOut = np.linspace(0, 100, 200)              # time
    # Initial case, 10 rabbits, 5 foxes
    Y0 = np.array([10, 5])
    t_span = [tOut[0], tOut[-1]] #tOut -1 stands gives back 100
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

#parameters landfill calculations 

baseA   = 28355     # unit m²
topA    = 9100      # unit m²
slopeW  = 38        # unit m
wastebH = 12        # unit m 
coverlH = 1.5       # unit m 
wasteWEI= 281083000 # unit kg
a       = 0.45      # unit m/d, saturated hydraulic conductivity, chosen from Shank (1993) 20year old landfill 
Sclmax  = coverlH * 0.35      # unit = m --> cover layer height * porosity, maximum acheivable storage in cover layer
#Swbmax  = XXX       # maximum achivable storage in waste body layer
Sclmin = 0.54      # unit m, minimum storage in cover layer where water will still freely drain
Swbmin  = 2.7       # unit m, minimum storage in waste body layer where water will still freely drain 
bcl     = 7.9       # dimensionless empirical parameter cover layer
bwb     = 28.0      # dimensionless empirical parameter waste body 

#test test change change
# define time dependent functions 

def dSdt(t, Y):
    """ Return the rates of change for the storage in the cover, waste and drainage layer"""
    return np.array([Jrf(rain, tspan) -  - Eva(evapo, tspan),
                                            ,  
                                            0])

def E(E, t):
    return E[t]

def J(J, t):
    return J[t]

def Lcl(Scl, t):
    a*(())



Eva(evapo, tspan)
Jrf(rain, tspan)

def mainLFWB:
    Eva = np.vectorize(E)
    Jrf  = np.vectorize(J)
    
    # read data from meteo file
    Meteo = pd.read_excel(r"C:\Users\tjark\OneDrive\Notability - Tjark\Modelling coupled processes CIE4365\Module 1 Intro & modelling Water Balance Landfill\WieringermeerData_Meteo.xlsx", index_col=0, parse_dates=[0])
    # Eva
    evapo = Meteo.pEV                 # unit m/d 
    #rain
    rain = Meteo.rain_station   # unit m/d
    
    
    # Definition of output times
    tOut = np.linspace(0, 6209, 6210)              # whole length of wetaher data 
    # Initial case, water stored in cover and drainage layer
    S0 = np.array([XX, XX])
    t_span = [tOut[0], tOut[-1]] #tOut -1 give whole length back 
    SODE = sp.integrate.solve_ivp(dSdt, t_span, S0, t_eval=tOut, 
                                  method='RK45', vectorized=True, 
                                  rtol=1e-5 )
    # infodict['message']                     # >>> 'Integration successful.'
    #rODE = YODE.y[0,:]
    #fODE = YODE.y[1,:]

    
    
Evaporation =   panda.read excel 
  