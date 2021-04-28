#!/usr/bin/env python
# coding: utf-8


import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import pandas as pd
#import MyTicToc as mt

# get_ipython().run_line_magic('matplotlib', 'inline')


def tic():
    #Homemade version of matlab tic and toc functions
    import time
    global startTime_for_tictoc
    startTime_for_tictoc = time.time()

def toc():
    import time
    if 'startTime_for_tictoc' in globals():
        print ("Elapsed time is " + str(time.time() - startTime_for_tictoc) + " seconds.")
    else:
        print ("Toc: start time not set")



# definition of parameters:
ba = 28355  # base area [m²] 
ta = 9100  # top area [m²]
sw = 38  # slope wdith [m]
wbh = 12  # waste body height [m]
clh = 1.5  # cover layer height [m]
w = 281083000  # waste (wet weight) [kg]

# Storage capacity = 1.5  # [mm/cm depth]
Scap = 1.5 * 0.1 * clh * ta  # [m³]
a = 1
fc = 0.6  # crop factor estimation, grassland
bc = 7.9  # from Benettin, dont know if right
bwb = 28  # from Benettin, dont know if right 
ß0 = 0.85  # from Benettin, dont know if right

'temporary assumed values'
Scmin = 0.1
Scmax = 0.9
Swbmin = 0.2
Swbmax = 0.9
Semin = 0.9
Semax = 0.5

M = pd.read_excel('WieringermeerData_Meteo.xlsx') # index_col='datetime'
# display(M)

pEv = M.pEV         # potential evapotranspiration [m/day] --> use this way: display(pEv['2003-01-03'])
R = M.rain_station  # inflow [m/day]
 


   
# Definition of rate equations:
'What values are predefined and what must be taken from the available data?'

def g(R, t):
    t = int(t)
    return R[t]
J = np.vectorize(g)

def f(pEv, S, t):
    if S[0] < Semin:
        fs = 0
    elif Semin <= S[0] <= Semax:
        fs = (S[0] - Semin) / (Semax - Semin)
    else: 
        fs = 1
    return pEv[t] * fc * fs
E = np.vectorize(f)

def h(S, t):
    return a * ((S[0] - Scmin) / (Scmax - Scmin)) ** bc
Lc = np.vectorize(h)

def dSdt(S, t):
    return J(R, t) - Lc(S, t) - E(pEv, S, t)

# def dSdt(t, S):
#    if S[0] < Semin:
#        fs = 0
#    elif Semin <= S[0] <= Semax:
#        fs = (S[0] - Semin) / (Semax - Semin)
#    else: 
#        fs = 1
#    Et = pEv[t] * fc * fs
#    Lc = a * ((S[0] - Scmin) / (Scmax - Scmin)) ** bc  # Sc[t]
 #   Lw = a * ((S[1] - Swbmin) / (Swbmax - Swbmin)) ** bwb  # Swb[t]
#    ß = ß0 * ((S[0] - Scmin) / (Scmax - Scmin))
#    Qdr = ß * Lc + Lw
#    J = Jt[t]

 #   return np.array([J - Lc - Et,
 #                   (1 - ß) * Lc - Lw])
 #                   # ß * Lc + Lw - Qdr])



def main():
# Definition of output times (0 to 2757 days)
    tOut = np.linspace(0, 1500, 2757)  
    # print(tOut)
    nOut = np.shape(tOut)[0]
    
    # Initial case: Sc = 0.7, Swb = 0.8
    S0 = np.array([0.7, 0.8])
    
    import scipy.integrate as spint
    
    tic()
    t_span = [tOut[0], tOut[-1]]
    YODE = spint.solve_ivp(dSdt, t_span, S0, t_eval=tOut, vectorized=True,
                           method='RK45', rtol=1e-5)
    # infodict['message']                     # >>> 'Integration successful.'
    ScODE = YODE.y[0,:]
    SwbODE = YODE.y[1,:]
    SdrODE = YODE.y[2,:]
    
    toc()
    
    plt.figure()
    plt.plot(tOut, ScODE, 'r-', label='Cover layer')
    plt.plot(tOut, SwbODE  , 'b-', label='Waste body')
    plt.plot(tOut, SdrODE  , 'g-', label='Drainage')
    
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.ylabel('Storage')
    plt.title('Storage rate landfill')

# f1.savefig('rabbits_and_foxes_1.png')

    plt.figure()
    plt.plot(SwbODE, ScODE, 'b-', label='ODE')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('waste body storage')
    plt.ylabel('Coverlayer storage')
    plt.title('Storage rate landfill')

if __name__ == "__main__":
    main()