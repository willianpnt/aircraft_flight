from functions import atmosisa, CAS2TAS, TAS2CAS, thrust
import numpy as np

class climb():
    def __init__(self, aircraft, h1, h2, xi=0, h_M=30000):
        '''
        INPUTS
        aircraft: Aircraft data
        h1: Start climb altitude (ft)
        h2: End climb altitude (ft)
        xi: Aircraft position when the climb starts (km)
        h_M: Altitude which the Mach number is known (ft)
        '''
        
        #--- Initilization initial data
        self.W = aircraft['W']            #(N)
        self.S = aircraft['S']            #(m2)
        self.T0 = aircraft['T0']          #(N)
        self.TSFC = aircraft['TSFC']      #(N/h/N)
        self.DISA = aircraft['DISA']      #(oC)
        self.CD0 = aircraft['CD0']        #(-)
        self.K = aircraft['K']            #(-)

        #---Initialization output data
        self.x = [xi]         #(km)
        self.time = 0         #(min)
        h_sclimb = [h2]       #(ft)

        #--- Compute calibrated airspeed
        # From 30000 ft, the climb is at Mach constant
        if h1 > 30000:              #Climb starting above 30000 ft -> mach constant
            self.M = aircraft['Mach']               #Mach number at starting climb (-)
        else:
            a = atmosisa(h_M, self.DISA)[1]
            self.M = aircraft['Mach']               #Mach number at altitude h_M (-)
            self.V = self.M*a                       #True airspeed at altitude h_M (m/s)                
            CAS = TAS2CAS(self.V, h_M, self.DISA)   #Calibrated airspeed at altitude h_M (= when the climb begins) (m/s)
        self.TAS = []

        #---Compute climb
        H = np.linspace(h1, h2, 10000)         #(ft)
        Hp = np.array([H[i]*atmosisa(H[i])[0]/atmosisa(H[i], self.DISA)[0]\
            for i in range(len(H))])         #Pressure altitude (ft)
        dH = Hp[1:]-Hp[:-1]                  #(ft)   

        for i in range(len(dH)):
            self.h = (H[i]+H[i+1])/2
            self.T = thrust(self.T0, self.h, self.DISA)       #Aircraft thrust (N)
            rho = atmosisa(self.h, self.DISA)[3]              #Air density (kg/m3)
            a = atmosisa(self.h, self.DISA)[1]                #Sound speed (m/s)
            if self.h <= 30000:                               #Below 30000 ft, climb at CAS constant
                self.V = CAS2TAS(CAS, self.h, self.DISA)      #True airspeed (m/s)
                self.M = self.V/a                             #Mach airspeed (-)
                facc = self.compute_facc(cte='CAS')           #Acceleration factor (CAS constant) (-)
            else:                                             #Above 30000 ft, climb at mach number constant
                self.V = self.M*a
                facc = self.compute_facc(cte='M')             #Acceleration factor (Mach cte) (-)
            self.TAS.append(self.V)                          
            RC = self.compute_RC(rho)/(1+facc)                #Rate of climb (fpm)
            if RC <= 500:                                     #Check stepclimb
                h_sclimb.append(self.h)
            dt = dH[i]/RC                                     #Time interval (min)
            self.time += dt                                   
            self.update_position(dt)                          #Update x position of the aircraft (km)
            self.update_W(dt)                                 #Update the aircraft weight (N)
        self.Wf = aircraft['W']-self.W                        #Weight of fuel consumed (N)
        self.h = (H[1:]+H[:-1])/2                            
        self.Range = self.x[-1] - self.x[0]                   #Range (km)
        self.x.pop(0)
        if len(h_sclimb) == 1: self.H_sc = h_sclimb[0]        #Stepclimb altitude (=h2 if don't need it) (ft)
        else: self.H_sc = h_sclimb[1]

    def compute_RC(self, rho):
        '''
        Compute the rate of climb
        '''
        q = 1/2*rho*self.V**2
        D = q*self.S*self.CD0 + self.K*self.W**2/(q*self.S)
        gamma = (self.T - D)/self.W
        RC = self.V*gamma*(60/0.3048)
        return RC
    
    def compute_facc(self, cte):
        '''
        Compute acceleration factor (= V/g*dV/dh)
        '''
        if self.h <= 11000/0.3048: 
            zeta = 0.190263*atmosisa(self.h)[0]/atmosisa(self.h, self.DISA)[0]
        else: zeta = 0
        if cte == 'CAS':
            psi = ((1+0.2*self.M**2)**3.5-1)/\
                (0.7*self.M**2*(1+0.2*self.M**2)**2.5) - zeta
        elif cte == 'M':
            psi = -zeta
        elif cte == 'EAS':
            psi = 1-zeta
        return 0.7*self.M**2*psi

    def update_W(self, dt):
        '''
        Update aircraft weight
        '''
        self.W -= self.TSFC*(dt*60)*self.T

    def update_position(self, dt):
        '''
        Update x position of aircraft
        '''
        self.x.append(self.x[-1] + self.V*(dt*60)/1000)

class cruise():
    def __init__(self, aircraft, h, par, par_type='Fuel', xi = 0):
        '''
        INPUTS
        aircraft: Aircraft data
        h: Cruize altitude (ft)
        par: Fuel percentage remaining at cruise's end (-) or cruise range (km)
        par_type =  If the parameter given (par) is fuel percentage or range
        xi: Aircraft position when the cruise starts (km)
        '''

        #--- Aircraft data
        DISA = aircraft['DISA']                           #(oC)
        rho = atmosisa(h, DISA, 'ft')[3]                  #(kg/m3)
        a = atmosisa(h, DISA, 'ft')[1]                    #(m/s)
        S = aircraft['S']                                 #(m2)
        CD0 = aircraft['CD0']                             #(-)
        K = aircraft['K']                                 #(-)
        self.V = aircraft['Mach']*a                       #(m/s)
        TSFC = aircraft['TSFC']                           #(N/h/N)
        self.W = aircraft['W']                            #(N)

        #--- Compute cruise
        #-- For fuel remaining percentage
        if par_type == 'Fuel':
            Wend = aircraft['BOW'] + par*aircraft['Wf_t'] #Fuel weight at cruise's end (N)
            W = np.linspace(self.W, Wend, 100)           
            self.Q = []                                   #Fuel flow
        
            for i in range(len(W)):
                CL = W[i]/(0.5*rho*self.V**2*S)           #Lift coefficient (-)
                CD = CD0 + K*CL**2                        #Drag coefficient (-)
                D = 1/2*rho*self.V**2*S*CD                #Drag force (N)
                self.Q.append(TSFC*D)
            
            r = np.array(self.V/self.Q)                   #Specific range (m/N)
            R = (r[1:]+r[:-1])/2*(W[:-1]-W[1:])                
            self.Range = np.sum(R)/1000                   #Cruise range (km)
            self.time = np.sum(R/self.V)/60               #Cruise time (min)
            self.W = Wend                                 
            self.Wf = aircraft['W'] - self.W              #Consumed fuel (N)
        
        #-- P/ alcance de cruzeiro
        elif par_type =='Distance':
            dt = 0.01                                     #Time interval (min)
            R = 0                                     
            t = 0                                     
            while (R < par):
                R += self.V*(dt*60)/1000                  #Travelled distance in dt (km)
                Cl = 2*self.W/(rho*S*self.V**2)
                Cd = CD0 + K*Cl**2
                T = 1/2*rho*self.V**2*S*Cd                #Engine thrust (N)
                self.W -= TSFC*(dt*60)*T                  #Update aircraft weight (km)
                t += dt
            self.Range = R                                #Cruise range (km)
            self.time = t                                 #Cruise time (min)
            self.Wf = aircraft['W'] - self.W              #Consumed fuel weight (N)
        self.x = [xi, xi + self.Range]
        self.h = [h, h]
        
class descent():
    def __init__(self, aircraft, h1, h2, xi = 0):
        '''
        INPUTS
        aircraft: Aircraft data
        h1: Altitude when the descent starts
        h2: Altitude when the descent ends
        xi: Aircraft position when the descent starts
        '''
        #--- Initilization initial data
        self.W = aircraft['W']            #(N)
        self.S = aircraft['S']            #(m2)
        self.T0 = aircraft['T0']          #(N)
        self.TSFC = aircraft['TSFC']      #(N/h/N)
        self.DISA = aircraft['DISA']      #(oC)
        self.CD0 = aircraft['CD0']        #(-)
        self.K = aircraft['K']            #(-)
        self.M = aircraft['Mach']         #(-)

        #--- Initialization output data
        self.x = [xi]            #Aircraft x position (min)
        self.time = 0            #(km)
        
        #--- Cálculo da velocidade calibrada
        a = atmosisa(h1, self.DISA)[1]            #Sound speed at h1 (m/s)
        self.V = self.M*a                         #True airspeed at beggining of the descent (m/s)
        CAS = TAS2CAS(self.V, h1, self.DISA)      #Clibrated airspeed at the descent (m/s)
        self.TAS = []

        #--- Compute descent
        H = np.linspace(h1, h2, 500)          #(ft)
        Hp = np.array([H[i]*atmosisa(H[i])[0]/atmosisa(H[i], self.DISA)[0]\
             for i in range(len(H))])         #Pressure altitude (ft)
        dH = Hp[1:]-Hp[:-1]                   #(ft)

        for i in range(len(dH)):
            self.h = (H[i]+H[i+1])/2
            self.T = thrust(self.T0, self.h, self.DISA)
            self.T *= 0.05                                #Thrust at the descent (5% of max) (N)
            rho = atmosisa(self.h, self.DISA)[3]          #Air density (kg/m3)
            a = atmosisa(self.h, self.DISA)[1]            #Sound speed (m/s)
            self.V = CAS2TAS(CAS, self.h, self.DISA)      #True airspeed (m/s)
            self.M = self.V/a                             #Mach number (-)
            facc = self.compute_facc('CAS')               #Acceleration factor (-)
            self.TAS.append(self.V)                       #True airspeed storage
            RD = self.compute_RD(rho)/(1+facc)            #Rate of sink (fpm)
            dt = dH[i]/RD                                 
            self.time += dt                               #Descent time (min)
            self.update_position(dt)                      #Range (km)
            self.update_W(dt)                             #Update aircraft weight (N)
        self.Wf = aircraft['W']-self.W                    #Consumed fuel weight (N)
        self.Range = self.x[-1] - self.x[0]
        self.x.pop(0)
        self.h = (H[1:]+H[:-1])/2                         #Descent altitudes (ft)
  
    def compute_RD(self, rho):
        '''
        Compute rate of sink
        '''
        q = 1/2*rho*self.V**2
        D = q*self.S*self.CD0 + self.K*self.W**2/(q*self.S)
        gamma = (self.T - D)/self.W
        RD = self.V*gamma*(60/0.3048)
        return RD
    
    def compute_facc(self, cte):
        '''
        Compute acceleration factor (= V/g*dV/dh)
        '''
        if self.h <= 11000/0.3048: 
            zeta = 0.190263*atmosisa(self.h)[0]/atmosisa(self.h, self.DISA)[0]
        else: zeta = 0
        if cte == 'CAS':
            psi = ((1+0.2*self.M**2)**3.5-1)/\
                (0.7*self.M**2*(1+0.2*self.M**2)**2.5)-zeta
        elif cte == 'M':
            psi = -zeta
        elif cte == 'EAS':
            psi = 1-zeta
        return 0.7*self.M**2*psi

    def update_W(self, dt):
        '''
        Update aircraft weight
        '''
        self.W -= self.TSFC*(dt*60)*self.T

    def update_position(self, dt):
        '''
        Update aircraft x position
        '''
        self.x.append(self.x[-1] + self.V*(dt*60)/1000)

class loiter():
    def __init__(self, aircraft, h, t_loiter, xi = 0):
        '''
        INPUTS
        aircraft: Aircraft data
        h: Loiter altitude (ft)
        t_loiter: Loiter time (min)
        xi: X position when the loiter start (km)
        '''
        
        #--- Initilization initial data
        self.W = aircraft['W']              #(N)
        S = aircraft['S']                   #(m2)
        DISA = aircraft['DISA']             #(oC)
        rho = atmosisa(h, DISA, 'ft')[3]    #(kg/m3)
        a = atmosisa(h, DISA, 'ft')[1]      #(m/s)
        TSFC = aircraft['TSFC']             #(N/h/N)
        CD0 = aircraft['CD0']               #(-)
        K = aircraft['K']                   #(-)
        self.V = aircraft['Mach']*a         #(m/s)
        
        #--- Compute loiter
        dt = 0.01                                 #Time interval (min)
        t = 0
        while t < t_loiter:
            Cl = 2*self.W/(rho*S*self.V**2)
            Cd = CD0 + K*Cl**2
            T = 1/2*rho*self.V**2*S*Cd            #Engine thrust at loiter (N)
            self.W -= TSFC*(dt*60)*T              #Update aircraft weight
            t += dt              
        self.time = t                             #Loiter time (min)
        self.Wf = aircraft['W'] - self.W          #Consumed fuel (N)
        self.x = xi                               #X position at the loiter (km)