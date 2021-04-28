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
import pandas as pd

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
    
    
# start land fill model
#parameters landfill calculations 
# take "educated guesses" as starting point 

P = 0.35            #soil porosity
baseA   = 28355     # unit m²
topA    = 9100      # unit m²
slopeW  = 38        # unit m
wastebH = 12        # unit m 
coverlH = 1.5       # unit m 
volCl   = coverlH * topA    # unit m³ cover layer
volWb   = wastebH * baseA   # unit m³ waste body  
wasteWEI= 281083000 # unit kg

Sclmax  = coverlH * P    # unit = m --> cover layer height * porosity, maximum acheivable storage in cover layer
Sclmin = Sclmax * 0.5       # unit m, minimum storage in cover layer where water will still freely drain
Swbmax  = wastebH * P    # maximum achivable storage in waste body layer
Swbmin  = Swbmax * 0.5      # unit m, minimum storage in waste body layer where water will still freely drain 
Cf = 0.8            # crop factor 


#parameterse to be optimized
a       = 0.45          # unit m/d, saturated hydraulic conductivity, chosen from Shank (1993) 20year old landfill 
bcl     = 7.9           # dimensionless empirical parameter cover layer
bwb     = 28.0          # dimensionless empirical parameter waste body 
beta0   = 0.85
SEmin   = 0.1 * Sclmax   # minimum storage of cl for evapo
SEmax   = 0.9 * Sclmax   # maximum storage of cl for evapo 

Meteo = pd.read_excel(r"C:\Users\tjark\OneDrive\Notability - Tjark\Modelling coupled processes CIE4365\Module 1 Intro & modelling Water Balance Landfill\WieringermeerData_Meteo.xlsx", index_col=0, parse_dates=[0])
# Eva
pev = np.array(Meteo.pEV)             # unit m/d 
#rain
rain = np.array(Meteo.rain_station)   # unit m/d

#Leachate production from measured values 
LPmeas = pd.read_excel(r"C:\Users\tjark\OneDrive\Notability - Tjark\Modelling coupled processes CIE4365\Module 1 Intro & modelling Water Balance Landfill\WieringermeerData_LeachateProduction.xlsx", index_col=0, parse_dates=[0])

# define the function to be solved

def dSdt(t, S):
    """ 
    Return the rates of change for the storage in the cover and waste layer
    
    S=[Scl, Swb]
    """
    
    # value correction before simulation; corrects values which are out of the defined range 
    if S[0] <  Sclmin:
        Scl = Sclmin
    elif S[0] > Sclmax:
        Scl = Sclmax
        
    if S[1] < Swbmin:
        Swb = Swbmin
    elif S[1] > Swbmax:
        Swb = Swbmax
        
    #reduction factor evapotranspiration
    if S[0] < SEmin:
        fred = 0
    elif S[0] > SEmax:
        fred = 1
    elif SEmin <= S[0] <= SEmax:
        fred = (S[0] - SEmin) / (SEmax - SEmin) 
        
    beta = beta0 * (S[0] - Sclmin) / (Sclmax - Sclmin) 
    Eva  = pev[int(t)] * Cf * fred
    Lcl  = a * ((S[0] - Sclmin) / (Sclmax - Sclmin)) ** bcl 
    Lwb  = a * ((S[1] - Swbmin) / (Swbmax - Swbmin)) ** bwb 
    
    
    return np.array([rain[int(t)] - Lcl - Eva, 
                     (1 - beta) * Lcl - Lwb])


def leach_correction_postcalc(SclODE, SwbODE):
    '''
    checks if calculated values are in the actual range of possible real values and corrects them if not 
    '''
    
    if SclODE[0] < Sclmin:
        SclODE = Sclmin
    elif SclODE[0] > Sclmax:
        SclODE = Sclmax
        
    if SwbODE[1] < Swbmin:
        SwbODE = Swbmin
    elif SwbODE[1] > Swbmax:
        SwbODE = Swbmax
    
    return SclODE, SwbODE
    
leach_correction_postcalc_vec = np.vectorize(leach_correction_postcalc)

def leachate_produc (SclODE, SwbODE):
    """
    this function calculates the quantity of the leachate production from the landfill Qdr in m³/d
    """
    
    beta = beta0 * (SclODE - Sclmin) / (Sclmax - Sclmin) 
    Lcl  = a * ((SclODE - Sclmin) / (Sclmax - Sclmin)) ** bcl 
    Lwb  = a * ((SwbODE - Swbmin) / (Swbmax - Swbmin)) ** bwb 
    Qdr  = beta * Lcl * topA + Lwb * baseA 
    return np.array(Qdr)
    

    


    
# Definition of output times
tOut = np.linspace(0, 6208, 6209)              # whole length of wetaher data 
# Initial case, water stored in cover and drainage layer
S0 = np.array([Sclmax * 0.2, Swbmax * 0.2])
t_span = [tOut[0], tOut[-1]] #tOut -1 give whole length back 
SODE = sp.integrate.solve_ivp(dSdt, t_span, S0, t_eval=tOut, 
                                  method='RK45', vectorized=True, 
                                  rtol=1e-5 )

SclODE = SODE.y[0, 3452:]
SwbODE = SODE.y[1, 3452:]
SclODE, SwbODE = leach_correction_postcalc_vec(SclODE, SwbODE)
 
# daily amount leachcate 
leachate = leachate_produc(SclODE, SwbODE)
# cumulative amount leachate 
cum_leachate = np.cumsum(leachate)
    
plt.clf()
plt.plot(LPmeas.index, LPmeas, label="measured values", color="red")
plt.plot(LPmeas.index, cum_leachate, label="model output", color="blue")
plt.legend()    
    
  