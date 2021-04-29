# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 17:33:02 2021

@author: Guido
"""

# matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
import scipy.optimize as opt


#Data import
import pandas as pd

meteo = pd.read_excel('WieringermeerData_Meteo.xlsx',
                      index_col=0, parse_dates=[0])
leachate = pd.read_excel('WieringermeerData_LeachateProduction.xlsx',
                         skiprows=0, parse_dates=[0],
                         names = ['Time', 'Leachate'])

# Definition of data
rain = np.array(meteo.rain_station)       #m/day
evap = np.array(meteo.pEV)                #m/day
temp = np.array(meteo.temp)               #celsius
leach_cum = np.array(leachate.Leachate)   #m3/day


# Definition of parameters
#Characteristics Cell VP-06 Wieringermeer landfill:
#waste base
A_wb = 28355        #[m2]
h_wb = 12           #[m]
V_wb = A_wb * h_wb  #[m3]
#cover layer
A_cl = 9100         #[m2]
h_cl = 1.5          #[m]
V_cl = A_cl * h_cl  #[m3]
#slope witdh [m]
w = 38
#waste (net weight [kg])
W = 281083000


p = 0.35                    #porosity 
S_cl_max = p * h_cl            
S_cl_min = 0.5 * S_cl_max           
S_wb_max = p * h_wb            
S_wb_min = 0.5 * S_wb_max         
C_f = 0.8
ESmin = 0.1 * S_cl_max
ESmax = 0.9 * S_cl_max

# a = 0.2
# b_cl = 7.9
# b_wb = 28
# beta_0 = 0.85
# S_EV_min = 0.0525
# S_EV_max = 0.4725


# Definition of Rate Equation
def water_balance(a, b_cl, b_wb, beta_0, S_EV_min, S_EV_max):
    def dYdt(t, Y):
        if Y[0] < S_cl_min:
            Y[0] = S_cl_min
        elif Y[0] > S_cl_max:
            Y[0] = S_cl_max
        
        if Y[1] < S_wb_min:
            Y[1] = S_wb_min
        elif Y[1] > S_wb_max:
            Y[1] = S_wb_max
  
    
        J_rf = rain
        f_red = 0
        if Y[0] < S_EV_min:
            f_red = 0
        elif Y[0] > S_EV_max:
            f_red = 1
        else:
            f_red = (Y[0] - S_EV_min) / (S_EV_max - S_EV_min)
        E = evap[int(t)] * C_f * f_red
        L_cl = a * ((Y[0] - S_cl_min) / (S_cl_max - S_cl_min))**b_cl
        L_wb = a * ((Y[1] - S_wb_min) / (S_wb_max - S_wb_min))**b_wb
        beta = beta_0 * ((Y[0] - S_cl_min) / (S_cl_max - S_cl_min))
        
        return np.array([J_rf[int(t)] - L_cl - E, (1 - beta) * L_cl - L_wb])
    
    
    def correction(S_cl, S_wb):
        if S_cl < S_cl_min:
            S_cl = S_cl_min
        elif S_cl > S_cl_max:
            S_cl = S_cl_max
        
        if S_wb < S_wb_min:
            S_wb = S_wb_min
        elif S_wb > S_wb_max:
            S_wb = S_wb_max
        
        return S_cl, S_wb
    
    correction_vec = np.vectorize(correction)
    
    def leachate(S_cl, S_wb):
        L_cl = a * ((S_cl - S_cl_min) / (S_cl_max - S_cl_min))**b_cl
        L_wb = a * ((S_wb - S_wb_min) / (S_wb_max - S_wb_min))**b_wb
        beta = beta_0 * ((S_cl - S_cl_min) / (S_cl_max - S_cl_min))
        
        return np.array(beta * L_cl * A_cl + L_wb * A_wb)

    tOut = np.linspace(0, 6208, 6209)
    t_span = [tOut[0], tOut[-1]]
    Y0 = np.array([S_cl_max * 0.1, S_wb_max * 0.1])
    
    YODE = integrate.solve_ivp(dYdt, t_span, Y0, t_eval=tOut,
                               method='RK45', vectorized=True,
                               rtol=1e-5)
    ScODE = YODE.y[0, 3452:]
    SwODE = YODE.y[1, 3452:]
    ScODE, SwODE = correction_vec(ScODE, SwODE)
    
    sim_data = leachate(ScODE, SwODE)
    cum_sim_data = np.cumsum(sim_data)
    
    plt.clf()
    plt.plot(tOut[:2757], leach_cum, label='Measured leachate')
    plt.plot(tOut[:2757], cum_sim_data, label='Simulated leachate')
    plt.title('Landfill leachate model fit')
    plt.legend()
    return 

water_balance(a = 0.2, b_cl = 8, b_wb = 28, beta_0 = 0.8, S_EV_min = 0.045, S_EV_max = 0.45)

    
        
    
    