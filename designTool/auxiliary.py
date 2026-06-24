# GENERAL IMPORTS
import numpy as np
from .constants import gamma_air, R_air

#=========================================

def atmosphere(z, deltaISA=0.0):

    '''
    Funçao que retorna a Temperatura, Pressao e Densidade para uma determinada
    altitude z [m]. Essa funçao usa o modelo padrao de atmosfera (ISA) para a
    temperatura no solo de 288.15+deltaISA.
    '''

    # Zbase (so para referencia)
    # 0 11019.1 20063.1 32161.9 47350.1 50396.4

    # DEFINING CONSTANTS
    # Earth radius
    r = 6356766
    # gravity
    g0 = 9.80665
    # air gas constant
    R = 287.05287
    # layer boundaries
    Ht = [0, 11000, 20000, 32000, 47000, 50000]
    # temperature slope in each layer
    A = [-6.5e-3, 0, 1e-3, 2.8e-3, 0]
    # pressure at the base of each layer
    pb = [101325, 22632, 5474.87, 868.014, 110.906]
    # temperature at the base of each layer
    Tstdb = [288.15, 216.65, 216.65, 228.65, 270.65]
    # air viscosity
    mi0 = 18.27e-6 # [Pa s]
    T0 = 291.15 # [K]
    C = 120 # [K]

    # geopotential altitude
    H = r*z/(r+z)

    # selecting layer
    if H < Ht[0]:
        raise ValueError('Under sealevel')
    elif H <= Ht[1]:
        i = 0
    elif H <= Ht[2]:
        i = 1
    elif H <= Ht[3]:
        i = 2
    elif H <= Ht[4]:
        i = 3
    elif H <= Ht[5]:
        i = 4
    else:
        raise ValueError('Altitude beyond model boundaries')

    # Calculating temperature
    T = Tstdb[i]+A[i]*(H-Ht[i])+deltaISA

    # Calculating pressure
    if A[i] == 0:
        p = pb[i]*np.exp(-g0*(H-Ht[i])/R/(Tstdb[i]+deltaISA))
    else:
        p = pb[i]*(T/(Tstdb[i]+deltaISA))**(-g0/A[i]/R)

    # Calculating density
    rho = p/R/T

    # Calculating viscosity with Sutherland's Formula
    mi=mi0*(T0+C)/(T+C)*(T/T0)**(1.5)

    # Calculating the speed of sound
    a = np.sqrt(gamma_air*R_air*T)

    # Store results in a dictionary
    atm_data = {'altitude':z,
                'deltaISA':deltaISA,
                'speed_of_sound':a,
                'temperature':T,
                'pressure':p,
                'density':rho,
                'dyn_viscosity':mi}

    return atm_data

#----------------------------------------

def lin_interp(x0, x1, y0, y1, x):
    '''
    Linear interpolation function
    '''

    y = y0 + (y1-y0)*(x-x0)/(x1-x0)

    return y