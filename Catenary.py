# -*- coding: utf-8 -*-
"""
Created on Thu Oct 27 17:59:28 2016

@author: Andrea
"""
import numpy as np
import matplotlib.pyplot as plt

class Catenary(object):
    def __init__(self):
        self.tension = 0
        self.payout = 0
        self.catRange = 0
        
    def setTension(self,tension):
        self.tension = tension
    def getTension(self):
        return self.tension
    def setPayout(self,payout):
        self.payout = payout
    def getPayout(self):
        return self.payout
    def setRange(self,catRange):
        self.catRange = catRange
    def getRange(self):
        return self.catRange
        
def calculateCatenary(tension,EA,weight,frld_height):
    """ 
    Calculate catenary given tension at fairlead, height of fairlead from
    seafloor and mooring line specification (weight per unit length and 
    stiffness EA.
        
    tension = tension of mooringg line at fairlead in [kN]
    EA = stiffness of the mooring line in [kN]
    weight = weight per unit length of the mooring line in [kN/m]
    frld_height = height of fairlead with respect to seabed in [m]
    
    x = list of longitudinal position of points in the catenary in [m]
    y = list of vertical position of points in the catenary in [m]
    theta = angle of the tangent to the catenary at specific points to the 
            horizontal measured in [rad]
    """    

### Formulae from 2.019 Design of Ocean Systems - Spring 2011 - MIT OpenCourseWare    
    t_horizontal = EA*np.sqrt((tension/EA+1)**2-(2*weight*frld_height/EA))-EA
    payout_susp = 1/weight*np.sqrt(tension**2-t_horizontal**2)
    payout_stretched = payout_susp*(1+tension/EA)
    t_vertical = weight*payout_susp                    

### Faltinsen formula - Page 270 "Sea Loads on Ships and Offshore Structures"                    
    min_range = t_horizontal/weight* \
        np.log(((t_horizontal**2+t_vertical**2)**0.5+t_vertical)/t_horizontal) \
                    +t_horizontal*payout_susp/EA
###
    theta0 = np.arccos(t_horizontal/tension) # angle at fairlead [radians]
###
### Formulae from GMoor Manual - Catenary Quasi-static theory 
    w_0 = weight*(1+t_horizontal/EA)
    cost = t_horizontal/w_0
    theta = []
    x = []
    y = []
    n_components = 21
    for i in range(n_components):
        theta.append(theta0/(n_components-1)*i)         
        x_theta = cost*np.log(1/np.cos(theta[i])+np.tan(theta[i]))+ \
                cost*np.tan(theta[i])*t_horizontal/EA
        x.append(x_theta)
        y_theta = cost*1/np.cos(theta[i])*\
            (1+0.5*t_horizontal/EA*1/np.cos(theta[i]))
        # Normalize y(s) with respect to the initial position (anchor) --> this is
        # done in order to include the constant of integration in the formula for
        # the evaluation of y(s)
        if i < 1:        
            y.append(y_theta)
            y0 = y[0]
            y[0] = y[0]-y0
        else:
            y.append((y_theta-y0))
            
    # Write output to screen to check results:
    print('Tension at fairlead [kN] = ', tension)
    print('Mooring line weight [kN/m] = ', weight)
    print('Mooring line stiffness EA [kN] = ',EA)
    print('Height of Fairlead [m] = ', frld_height,'\n')
    print('TH [kN] = ', t_horizontal)    
    print('TZ [kN] = ',t_vertical)
    print('Angle at fairlead [deg] = ' + str(theta0*180/np.pi))
    print('Payout suspendded [m] = ',payout_susp)
    print('Stretched payout [m] = ', payout_stretched)
    print('Range Fld-TDP [m] = ', min_range)   
#    print('X(theta) [m] = ', x_theta)
#    print('y[0] [m] = ', y[0])
#    print('y[-1] [m] = ', y[-1])
    
    plt.figure('Cat')
    plt.clf()
    plt.plot(x,y)
    plt.xlabel('X [m]')
    plt.ylabel('Y [m]')
    plt.title('Catenary - Top Tenstion = '+str(tension)+'kN')
    plt.grid(which='major')    
    
    return x, y, theta

""" Testing the functions """
# Initial Inputs:
EA = 10**9 # [kN]
weight = 0.828 # [kN]
frld_height = 100 # [m]
tension = 1510 # [kN]

# Test the function/method:
x, y, theta = calculateCatenary(tension,EA,weight,frld_height)

# Print Results:
for i in range(len(x)):
    print(x[i],y[i],theta[i])