import numpy as np
import matplotlib.pyplot as plt

def atmosisa(h, D_ISA = 0, unit='ft'):
    '''
    Compute international standard atmosphere

    INPUT:
    h: Altitude (ft)
    DISA: Delta ISA
    unit: Altitude unit (ft ou m) 
   '''
    if unit == 'ft':
        h *= 0.3048
    c = -6.5e-3
    R = 287.053
    gamma = 1.4
    T_ssl = 288.15
    rho_ssl = 1.225
    g =  9.80665
    if h <= 11000:
        T = T_ssl + c*h
        rho = rho_ssl*(T/T_ssl)**(-g/(c*R)-1)
    else:
        T = T_ssl + c*11000
        rho = rho_ssl*(T/T_ssl)**(-g/(c*R)-1)\
            *np.exp(-g*(h-11000)/(R*T))
    p = rho*R*T
    T += D_ISA
    rho = p/(R*T)
    a = np.sqrt(gamma*R*T)
    return T, a, p, rho

def thrust(T0, h, DISA = 0):
    '''
    Compute the engine thrust at a given altitude
    
    INPUT:
    T0: Thrust at sealevel 
    h: Altitude
    DISA: Delta isa
    '''
    rho0 = atmosisa(0)[3]
    rho = atmosisa(h, DISA)[3]
    sigma = rho/rho0
    if h <= 11000/0.3048:
        T = T0*sigma**0.7
    else:
        T = 1.439*T0*sigma
    return T

def Mach_long_range(aircraft, h, DISA = 0):
    '''
    Compute the mach number for long range cruise

    INPUT:
    aircraft: Aircraft data
    h: Altitude
    DISA: Delta ISA
    '''
    V_vec = np.linspace(100, 300, 1000) 
    T = thrust(aircraft['T0'], h, DISA)
    a = atmosisa(h, DISA, 'ft')[1] 
    rho = atmosisa(h, DISA, 'ft')[3]
    M_vec = V_vec/a
    q = 1/2*rho*V_vec**2
    T = q*aircraft['S']*aircraft['CD0'] + aircraft['K']*aircraft['W']**2/(q*aircraft['S'])
    SR = V_vec/(aircraft['TSFC']*T)
    #--- Plot specific range
    plt.figure()
    plt.plot(M_vec, SR*1.852)
    plt.grid()
    plt.title('Specific range')
    plt.xlabel('M [-]')
    plt.ylabel('SR [nm/ton]')
    plt.savefig('SRMach.png')
    
    p = np.where(SR == max(SR))[0][0]
    SR_LR = 0.99*max(SR)
    M = np.interp(SR_LR, SR[:p+1] ,M_vec[:p+1]) 
    return M

def CAS2TAS(CAS, h, DISA):
    '''
    Convert calibrated airspeed in true airspeed

    INPUT:
    CAS: Calibrated airspeed(m/s)
    h: Altitude (ft)
    DISA: Delta ISA
    '''
    T, a, P, rho = atmosisa(h, DISA)
    T0, a0, P0, rho0 = atmosisa(0)
    EAS = np.sqrt(2.8/0.4*(P/rho0))\
        *np.sqrt((((1+(0.4/2)*(CAS/a0)**2)**(1.4/0.4)-1)/(P/P0)+1)**(0.4/1.4)-1)
    TAS = EAS/np.sqrt(rho/rho0)
    return TAS

def TAS2CAS(TAS, h, DISA):
    '''
    Convert true airspeed in calibrated airspeed

    INPUT:
    TAS: True airspeed (m/s)
    h: Altitude (ft)
    DISA: Delta ISA
    '''
    T, a, P, rho = atmosisa(h, DISA)
    T0, a0, P0, rho0 = atmosisa(0)
    EAS = TAS*np.sqrt(rho/rho0)
    CAS = a0*np.sqrt((((((EAS**2)/(2.8/0.4*(P/rho0))+1)**(1.4/0.4)-1)*(P/P0)+1)**(0.4/1.4)-1)*(2/0.4))
    return CAS