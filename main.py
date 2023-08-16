import numpy as np
import matplotlib.pyplot as plt
import flight as f   
from functions import atmosisa, Mach_long_range

#%% Aircraft data

aircraft = {}
aircraft['DISA'] = +15
aircraft['AR'] = 8.9
aircraft['S'] = 92.5
aircraft['e'] = 0.85
aircraft['CD0'] = 0.025
aircraft['K'] = 1/(aircraft['e']*aircraft['AR']*np.pi)
aircraft['MTOW'] = 450300
aircraft['BOW'] = 320300
aircraft['Wf_t'] = 130000
aircraft['T0'] = 92300
aircraft['TSFC'] = 0.85/3600
aircraft['Wf'] = aircraft['Wf_t']
aircraft['W'] = aircraft['MTOW']
aircraft['Mach'] = 0.79

#-- Plots
def plot_h_x(mission):
    try: plt.plot(mission.x, mission.h, color='k')
    except: pass

def plot_h_V(mission):
    try: plt.plot(mission.TAS, mission.h, color='k')
    except: pass
#%% Flight info

FL = 41000                             #Cruise flight level (ft)
FL_alt = 20000                         #Alternative cruise flight level (ft)
step = 2000
FL_loiter = 10000                      #Loiter flight level (ft)
Loiter_time = 30                       #Loiter time (min)

h_airport_dp = 0                       #Departure airport altitude (ft)
h_airport_alt = 2000                   #Alternative airport altitude (ft)
h_airport_dest = 0                     #Arrival airport altitude (ft)

#%% Initialization 

Wf_c = 0
mission = []

#%% Climb
print('*** Climb phase ---------')

mission.append(f.climb(aircraft, h_airport_dp, FL))
h_sc = FL - step*np.ceil((FL - mission[-1].H_sc)/step)

if h_sc != FL:
    print('The aircraft reaches R/C = 500 fpm at %d ft' %mission[-1].H_sc)
    print('Need a stepclimb at %d ft\n' %h_sc)
    x_sc = 100
    h0 = [0]
    mission.pop()
    while h_sc != FL:
        try: xi = mission[-1].x[-1]
        except: xi = 0
        W_0 = aircraft['W']
        mission.append(f.climb(aircraft,h0[-1],h_sc, xi))
        aircraft['W'] = mission[-1].W
        aircraft['Mach'] = mission[-1].M
        mission.append(f.cruise(aircraft, h_sc, x_sc, 'Distance',mission[-1].x[-1]))
        aircraft['W'] = mission[-1].W
        mission.append(f.climb(aircraft, h_sc, FL, mission[-1].x[-1]))
        h0.append(h_sc)
        h_sc = FL - step*np.ceil((FL - mission[-1].H_sc)/step)
        if h0[-1] == h_sc:
            x_sc += 100
            h0.pop()
            mission.pop()
            mission.pop()
            mission.pop()
            aircraft['W'] = W_0
        else:
            Wf_c += mission[-3].Wf
            Wf_c += mission[-2].Wf
            print('-> Climb until FL%d ----' %(h_sc/100 - step/100))
            print('\tTime to climb: %.2f min' %(mission[-3].time))
            print('\tTravelled distance: %.2f km' %mission[-3].Range)
            print('\tConsumed fuel: %.2f N\n' %mission[-3].Wf)
            print('-> Stepclimb at FL%d ---' %(h_sc/100 - step/100))
            print('\tTime of stepclimb:  %.1f min' %mission[-2].time)
            print('\tTravelled distance: %.2f km' %mission[-2].Range)
            print('\tConsumed fuel: %.2f N\n' %mission[-2].Wf)
            if h_sc != FL:
                mission.pop()
    aircraft['W'] = mission[-1].W
    aircraft['Wf'] = aircraft['W'] - aircraft['BOW']
    Wf_c += mission[-1].Wf
    print('-> Climb until FL%d ----' %(h_sc/100))
    print('\tTime to climb: %.2f min' %(mission[-1].time))
    print('\tTravelled distance: %.2f km' %mission[-1].Range)
    print('\tConsumed fuel: %.2f N\n' %mission[-1].Wf)
else:
    Wf_c += mission[-1].Wf
    print('-> Dont need a stepclimb <-')
    print('\tTime to climb: %.1f min' %mission[-1].time)
    print('\tTravelled distance: %.1f km' %mission[-1].Range)
    print('\tConsumed fuel: %.2f N\n' %mission[-1].Wf)
    aircraft['W'] = mission[-1].W
    aircraft['Wf'] = aircraft['W'] - aircraft['BOW']

plt.figure()
for i in range(len(mission)):
    plot_h_V(mission[i])
plt.xlabel('TAS [m/s]')
plt.ylabel('h [ft]')
plt.grid()
    
#%% Cruise
print('*** Cruise phase --------')

mission.append(f.cruise(aircraft, FL, 0.19, xi= mission[-1].x[-1]))
print('\tTime in cruise: %.1f min' %mission[-1].time)
print('\tTravelled distance: %.1f km' %mission[-1].Range)
print('\tConsumed fuel: %.2f N\n' %mission[-1].Wf)

Wf_c += mission[-1].Wf
aircraft['W'] = mission[-1].W
aircraft['Wf'] = aircraft['W'] - aircraft['BOW']


#%% Descent
print('*** Descent phase -------')

mission.append(f.descent(aircraft, FL, h_airport_dest, mission[-1].x[-1]))
print('\tTime to descent: %.1f min' %mission[-1].time)
print('\tTravelled distance: %.1f km' %mission[-1].Range)
print('\tConsumed fuel: %.2f N\n' %mission[-1].Wf)

plt.figure()
plt.plot(mission[-1].TAS, mission[-1].h)
plt.xlabel('TAS [m/s]')
plt.ylabel('h [ft]')
plt.grid()

Wf_c += mission[-1].Wf
Wf_dest = Wf_c
aircraft['W'] = mission[-1].W
aircraft['Wf'] = aircraft['W'] - aircraft['BOW']

#-- Aircraft range
R = 0                                           #Aircraft range (km)
E = 0                                           #Aicraft endurance (min)
for i in range(len(mission)):
    E += mission[i].time
    R += mission[i].Range

#%% Climb to alternative
print('*** Climb to alternative ---------')

x_alt = 0
aircraft['DISA'] = +20
aircraft['Mach'] = Mach_long_range(aircraft, FL_alt, aircraft['DISA'])
Cross_h = 20000

mission.append(f.climb(aircraft, h_airport_dest, FL_alt, mission[-1].x[-1], Cross_h))

Wf_c += mission[-1].Wf
x_alt += mission[-1].Range
print('\tTime to climb: %.1f min' %mission[-1].time)
print('\tTravelled distance: %.1f km' %mission[-1].Range)
print('\tConsumed fuel: %.2f N\n' %mission[-1].Wf)
aircraft['W'] = mission[-1].W
aircraft['Wf'] = aircraft['W'] - aircraft['BOW']
#%% Cruise

print('*** Cruise in alternative phase --------')

mission.append(f.cruise(aircraft, FL_alt, 230, 'Distance', mission[-1].x[-1]))
print('\tTime in cruise: %.1f min' %mission[-1].time)
print('\tTravelled distance: %.1f km' %mission[-1].Range)
print('\tConsumed fuel: %.2f N\n' %mission[-1].Wf)

Wf_c += mission[-1].Wf
x_alt += mission[-1].Range
aircraft['W'] = mission[-1].W
aircraft['Wf'] = aircraft['W'] - aircraft['BOW']

#%% Descent 1
print('*** Descent to loiter altitude -------')

mission.append(f.descent(aircraft, FL_alt, FL_loiter, mission[-1].x[-1]))
print('\tTime to descent: %.1f min' %mission[-1].time)
print('\tTravelled distance: %.1f km' %mission[-1].Range)
print('\tConsumed fuel: %.2f N\n' %mission[-1].Wf)

Wf_c += mission[-1].Wf
x_alt += mission[-1].Range
aircraft['W'] = mission[-1].W
aircraft['Wf'] = aircraft['W'] - aircraft['BOW']

#%% Loiter

print('*** Loiter phase ----------------------')

rho = atmosisa(FL_loiter, aircraft['DISA'], 'ft')[3]
a = atmosisa(FL_loiter, aircraft['DISA'], 'ft')[1]
V = np.sqrt(2/rho*aircraft['W']/aircraft['S']\
            *np.sqrt(aircraft['K']/aircraft['CD0']))
aircraft['Mach'] = V/a

mission.append(f.loiter(aircraft, FL_loiter, 30, 'Time'))
print('\tTime in loiter: %.1f min' %mission[-1].time)
print('\tConsumed fuel: %.2f N\n' %mission[-1].Wf)

Wf_c += mission[-1].Wf
aircraft['W'] = mission[-1].W
aircraft['Wf'] = aircraft['W'] - aircraft['BOW']

#%% Descent 2

print('*** Descent to alternative airport -------')

mission.append(f.descent(aircraft, FL_loiter, h_airport_alt, mission[-2].x[-1]))
print('\tTime to descent: %.1f min' %mission[-1].time)
print('\tTravelled distance: %.1f km' %mission[-1].Range)
print('\tConsumed fuel: %.2f N\n' %mission[-1].Wf)

Wf_c += mission[-1].Wf
x_alt += mission[-1].Range
aircraft['W'] = mission[-1].W
aircraft['Wf'] = aircraft['W'] - aircraft['BOW']

#%% Plots

plt.figure()
for i in range(len(mission)):
    plot_h_x(mission[i])
plt.xlabel('x [km]')
plt.ylabel('h [ft]')
plt.title('Aircraft flight mission')
plt.grid()
print('missao.png')

print('Travelled distance until alternative airport: %.2f nm' %(x_alt*0.54))
print('Aircraft range : %.1f km' %R)
print('Flight time: %.1f h' %(E/60))
print('Remaining fuel: %.2f N' %(aircraft['Wf']))
print('Reserve fuel: %.2f' %((Wf_c - Wf_dest)*100/aircraft['Wf_t']))