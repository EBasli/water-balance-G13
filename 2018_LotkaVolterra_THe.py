#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import pandas as pd
import scipy.integrate as spint
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
ba = 28355              # base area [m²] 
ta = 9100               # top area [m²]
sw = 38                 # slope wdith [m]
wbh = 12                # waste body height [m]
clh = 1.5               # cover layer height [m]
w = 281083000           # waste (wet weight) [kg]

p = 0.35                # from Benettin
fc = 0.6                # crop factor estimation, grassland

'where does the 0.5 come from - what source?'
Scmin = ta * clh * p * 0.5  
Scmax = ta * clh * p
Swbmin = ba * wbh * p * 0.5
Swbmax = ba * wbh * p 

'value estimations - to be adjusted'
a = 741                # starting value 0.48 from Benettin
bc = 0.26                # 7.9 from Benettin
bwb = 100                # 28 from Benettin
ß0 = 26.5              # 0.85 from Benettin
Semin = 0.2 * Scmax     # start 0.1 * Scmax
Semax = 0.9 * Scmax     # start 0.9 * Scmax

M = pd.read_excel('WieringermeerData_Meteo.xlsx') # index_col='datetime'
L = pd.read_excel('WieringermeerData_LeachateProduction.xlsx')

pEv = M.pEV                     # potential evapotranspiration [m/day] --> use this way: display(pEv['2003-01-03'])
J = M.rain_station              # inflow [m/day]
# headers edited in excel file
LP = L.Leachate_Production      # Leachate Production [m³/day]

   
# Definition of rate equations:
def dSdt(t, S):
    # making sure storage value is within boundries:
    if S[0] < Scmin:
        S[0] = Scmin
    elif S[0] > Scmax:
        S[0] = Scmax
    if S[1] < Swbmin:
        S[1] = Swbmin
    elif S[1] > Swbmax:
        S[1] = Swbmax 
    
    # rate equations:
    if S[0] < Semin:
        fr = 0
    elif Semin <= S[0] <= Semax:
        fr = (S[0] - Semin) / (Semax - Semin)
    else: 
        fr = 1
    Et = pEv[int(t)] * fc * fr
    Lc = a * ((S[0] - Scmin) / (Scmax - Scmin)) ** bc  # Sc[t]
    Lw = a * ((S[1] - Swbmin) / (Swbmax - Swbmin)) ** bwb  # Swb[t]
    ß = ß0 * ((S[0] - Scmin) / (Scmax - Scmin))

    dSc = J[int(t)] - Lc - Et
    dSwb = (1 - ß) * Lc - Lw
    return np.array([dSc, dSwb])
                   # ß * Lc + Lw - Qdr])

def value_correction(ScODE, SwbODE):
    if ScODE < Scmin:
        ScODE = Scmin
    elif ScODE > Scmax:
        ScODE = Scmax
    if SwbODE < Swbmin:
        SwbODE = Swbmin
    elif SwbODE > Swbmax:
        SwbODE = Swbmax
    
    Sc_new = ScODE
    Swb_new = SwbODE
    return Sc_new, Swb_new

value_correction_vectorize = np.vectorize(value_correction)

def L_production(Sc_new, Swb_new):
    Lc = a * ((Sc_new - Scmin) / (Scmax - Scmin)) ** bc * ta # Sc[t]
    Lw = a * ((Swb_new - Swbmin) / (Swbmax - Swbmin)) ** bwb * ba # Swb[t]
    ß = ß0 * ((Sc_new - Scmin) / (Scmax - Scmin))
    Qdr = ß * Lc + Lw
    return Qdr

def main():
# Definition of output times (0 to 2757 days)
    tOut = np.linspace(0, 6209, 6210)  
    #nOut = np.shape(tOut)[0]
    
    # Initial case: Sc = 0.7, Swb = 0.8
    S0 = np.array([Semin, Semax])
    
    
    
    tic()
    t_span = [tOut[0], tOut[-1]]
    SODE = spint.solve_ivp(dSdt, t_span, S0, t_eval=tOut, vectorized=True,
                           method='RK45', rtol=1e-5)
    # infodict['message']                     # >>> 'Integration successful.'
    ScODE = SODE.y[0, 3452:]
    SwbODE = SODE.y[1, 3452:]
    Sc_new, Swb_new = value_correction_vectorize(ScODE, SwbODE)
    
    Qdr_day = L_production(Sc_new, Swb_new) # total daily leachate production
    Qdr_tot = np.cumsum(Qdr_day) # total leachate production in time period defined
    #QF = np.sum(Qdr_tot - LP) ** 2
    
    
    
    toc()
    
    plt.figure()
    plt.plot(LP.index, LP, 'r-', label='leachate measured')
    plt.plot(range(len(Qdr_tot)), Qdr_tot  , 'b-', label='leachate simulated total')
    #plt.plot(tOut, Qdr_day, 'g-', label='leachate simulated daily')
    #plt.plot(tOut, Sc_new, 'l-', label='storage cover layer')
    #plt.plot(tOut, Swb_new, 'p-', label='storage waste body')
    
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('time')
    plt.ylabel('Storage')
    plt.title('Storage rate landfill')

# f1.savefig('rabbits_and_foxes_1.png')

    plt.figure()
    plt.plot(Swb_new, Sc_new, 'bo', label='ODE')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('waste body storage')
    plt.ylabel('Coverlayer storage')
    plt.title('Storage rate landfill')
    
    #return QF

if __name__ == "__main__":
    main()