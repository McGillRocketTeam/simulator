# -*- coding: utf-8 -*-
"""
McGill Rocket Team Flight Simulator
Authors: Matt Saathoff, Mei Qi Tang

Internally, this simulator uses only SI units. When inputting values, it is currently up to the user to make sure
the units are correct. At some point we should make a GUI that supports unit conversion when inputting values,
but we should still do all our calculations in SI.

There is still a lot of structure left unmade too, as well as variables which need to be passed between functions.
I can keep working on that, but at some point I'm not totally sure what each function will need. Also, each
variable that is passed to a creator function needs to be checked for improper inputs to avoid running into
python's native exception handling. I've been doing this with if statements, print statements, and error state
attributes, but if someone wants convert this to proper try/catch exception statements that would probably be good.
"""

#Import Libraries
import numpy as np
import matplotlib.pyplot as plt
import IPython.display as IPy
from math import factorial
import helper

#A Class for Creating Rocket Objects
class Rocket:
    #A Class for Creating Motor Objects
    class Motor:
        #List of Supported Motor Types
        motorTypes = ['commercial solid','solid','hybrid','liquid']
        
        #Usefull Constants
        R = 8.3144598   #J/K*mol
        
        #A Class for Creating Tank Objects
        class Tank:
            #Lists of Supported Propellants
            liquidFuels = ['RP1','ethanol','methanol','LNG','gasoline','N2H4','MMH','UDMH']
            hybridFuels = ['paraffin','HTPB','PE']
            oxidizers = ['N2O','LOx','nitrox','H2O2','ClF3','ClF5','N2O4','N2F4','ClO3F','FClO4','FLOx']
            pressurants = ['He','N2','Ar']
            
            #Create a tank Object
            def __init__(self, motorType, propellant, volume, propellantMass, initialPressure, outletDiameter, pressurant = None, LOxFraction = None, verbose = False):
                #Error State Attribute
                self.error = False
                
                #Check/Assign the Propellant Type/Key
                if motorType == 'liquid' and propellant in self.liquidFuels:
                    self.type = 'liquid fuel'
                    self.propellant = propellant
                    for key in range(len(self.liquidFuels)):
                        if propellant == self.liquidFuels[key]:
                            self.key = key
                            break
                elif motorType == 'hybrid' and propellant in self.hybridFuels:
                    self.type = 'hybrid fuel'
                    self.propellant = propellant
                    for key in range(len(self.hybridFuels)):
                        if propellant == self.hybridFuels[key]:
                            self.key = key
                            break
                elif motorType in {'liquid','hybrid'} and propellant in np.append(self.liquidFuels, self.hybridFuels):
                    print('Invalid Fuel Selection')
                    if motorType == 'liquid':
                        print('Supported Liquid Fuels: ' + str(self.liquidFuels))
                    else:
                        print('Supported Hybrid Fuels: ' + str(self.hybridFuels))
                    self.error = True
                    return
                elif propellant in self.oxidizers:
                    self.type = 'oxidizer'
                    self.propellant = propellant
                    for key in range(len(self.oxidizers)):
                        if propellant == self.oxidizers[key]:
                            self.key = key
                            break
                    if propellant in {'nitrox','FLOx'}:
                        if type(LOxFraction) not in {float, np.float64} or LOxFraction < 0 or LOxFraction > 1:
                            print('LOxFraction must be a number between 0 and 1')
                            self.error = True
                            return
                else:
                    print('Propellants Not Supported')
                    if motorType == 'liquid':
                        print('Supported Liquid Fuels: ' + str(self.liquidFuels))
                    else:
                        print('Supported Hybrid Fuels: ' + str(self.hybridFuels))
                    print('Supported Oxidizers: ' + str(self.oxidizers))
                    self.error = True
                    return
                
                #Check the Volume
                if helper.validateNumber(volume):
                    self.V = volume
                else:
                    print('Volume must be a number greater than 0')
                    self.error = True
                    return

                #Check/Assign the Propellant Mass
                if helper.validateNumber(propellantMass):
                    self.m = np.array([propellantMass])
                else:
                    print('Propellent Mass must be a number greater than 0')
                    self.error = True
                    return
                
                #Check/Assign the Initial Pressure
                if helper.validateNumber(initialPressure):
                    self.P = np.array([initialPressure])
                else:
                    print('Initial Pressure must be a number greater than 0')
                    self.error = True
                    return
                
                #Check/Assign the Outlet Diameter
                if helper.validateNumber(outletDiameter):
                    self.D = outletDiameter
                else:
                    print('Outlet Diameter must be a number greater than 0')
                    self.error = True
                    return
                
                #Check/Assign the Pressurant
                if pressurant and pressurant in self.pressurants:
                    self.pressurant = pressurant
                else:
                    if type(pressurant) == str:
                        print('Pressurant Not Supported')
                        print('Supported Pressurants: ' + str(self.pressurants))
                        self.error = True
                        return
                    elif verbose == True:
                        print('Invalid Pressurant')
                    self.error = True
                
                #Assign the Verbose Attribute
                self.verbose = verbose
            
            #Initialize a Tank Before Launch
            def initialize(self, temperature):
                pass
            
            #Update a Tank
            def updateTank(self, a, theta):
                pass
                
        #A Class for Creating Combustion Chamber Objects
        class Chamber:
            #Create a chamber Object
            def __init__(self, motorType, ID, length, fuelLength = None, portDiameter = None, preCombustion = None, verbose = False):
                #Error State Attribute
                self.error = False
                
            
            #Solid Motor Combustion Chamber Model
            def updateSolid(self, h):
                return None, None
            
            #Hybrid Motor Combustion Chamber Model
            def updateHybrid(self, ODot, h):
                return None, None
            
            #Liquid Motor Combustion Chamber Model
            def updateLiquid(self, ODot, FDot, h):
                return None, None
        
        #A Class for Creating Nozzle Objects
        class Nozzle:
            #List of Supported Nozzle Types
            nozzleTypes = ['CD']
            
            #Create a Nozzle Object
            def __init__(self, nozzleType, inletRadius, throatRadius, areaRatio, outletAngle, k, molecularWeight, efficiency = 1, verbose = False):
                #Check the Nozzle Type and get its Key
                if nozzleType in self.nozzleTypes:
                    #Find the Nozzle Key
                    for key in range(len(self.nozzleTypes)):
                        if self.nozzleTypes[key] == nozzleType:
                            self.key = key
                            break
                else:
                    print('Nozzle Type Not Supported')
                    print('Supported Nozzle Types: ' + str(self.nozzleTypes))
                    self.error = True
                    return
                
                #Check the Throat Radius
                if not helper.validateNumber(throatRadius):
                    print('Throat Radius must be a number greater than 0')
                    self.error = True
                    return
                
                #Check the Inlet Radius
                if not helper.validateNumber(inletRadius) or inletRadius <= throatRadius:
                    print('Inlet Radius must be a number greater than the Throat Radius')
                    self.error = True
                    return
                
                #Check the Area Ratio
                if not helper.validateNumber(areaRatio) or areaRatio <= 1:
                    print('Area Ratio must be a number greater than 1')
                    self.error = True
                    return
                
                #Check the Outlet Angle
                if not helper.validateNumber(outletAngle) or outletAngle < 0 or outletAngle > 90:
                    print('Outlet Angle must be a number between 0 and 90 degrees')
                    self.error = True
                    return
                
                #Check k
                if not helper.validateNumber(k) or k < 1:
                    print('k must be a number greater than 1')
                    self.error = True
                    return
                
                #Check the Molecular Weight
                if not helper.validateNumber(molecularWeight) or molecularWeight < 0:
                    print('Molecular weight must be a number greater than 0')
                    self.error = True
                    return
                
                #Check the Efficiancy
                if not helper.validateNumber(efficiency) or efficiency > 1:
                    print('Efficiency must be a number between 0 and 1')
                    self.error = True
                    return
                
                #Assign All Attributes
                self.nozzleType = nozzleType
                self.epsilon = areaRatio
                self.ri = inletRadius
                self.rt = throatRadius
                self.re = throatRadius*np.sqrt(areaRatio)
                self.theta = outletAngle
                self.k = k
                self.M = molecularWeight
                self.R = Rocket.Motor.R
                self.F = 0
                self.mDot = 0
                self.efficiency = efficiency
                self.choked = False
                self.error = False
                self.verbose = verbose
            
            #Update the Nozzle (Assuming Isentropic Flow)
            def updateNozzle(self, P0, T0, Patm, precision = 10**-10):
                #Check for Choked Flow
                self.choked = (Patm <= P0*(2/(self.k + 1))**((self.k + 1)/(self.k - 1)))
                
                #Determine the Thrust and Mass Flow
                if self.choked:
                    #Calculate the Mass FLow through the Throat
                    self.mDot = self.efficiency*P0*np.pi*(self.rt**2)*np.sqrt(self.k/(self.R*T0))*(1 + (self.k - 1)/2)**((self.k + 1)/(2 - 2*self.k))
                    
                    #Calculate the Exit Mach
                    """
                    There is almost certainly a better way to solve this, but this method isn't bad
                    either since it runs in O(log(Me/precision))
                    """
                    Max = 1
                    while self.epsilon > ((1 + 0.5*(self.k - 1)*Max**2)**(0.5*(self.k + 1)/(self.k - 1)))/Max:
                        Max *= 2
                    Min = Max/2
                    while abs(Max - Min) > precision:
                        Me = (Min + Max)/2
                        if self.epsilon > ((1 + 0.5*(self.k - 1)*Me**2)**(0.5*(self.k + 1)/(self.k - 1)))/Me:
                            Min = Me
                        else:
                            Max = Me

                    Me = (Min + Max) / 2 #not sure about this here...
                    #Calculate Pe/P0 (Exhaust Pressure over Stagnation Pressure)
                    Pr = ((2/(self.k + 1))/(1 + 0.5*(self.k - 1)*Me**2))**(self.k/(self.k - 1))
                    
                    #Calculate the Exhaust Velocity
                    Ve = np.sqrt(2*self.R*T0*self.k*((1 - (Pr)**((self.k - 1))/self.k))/(self.M*(self.k - 1)))
                    
                    #Calculate the Thrust
                    self.F = self.mDot*Ve + (P0*Pr - Patm)*np.pi*self.re**2
                else:
                    pass
        
        #Create a Motor Object
        def __init__(self, motorType, fuel = None, fuelMass = None, fuelVolume = None, fuelPressure = None, fuelOutletDiameter = None, fuelPressurant = None, oxidizer = None, oxMass = None, oxVolume = None, oxPressure = None, oxOutletDiameter = None, oxPressurant = None, LOxFraction = None, chamberID = None, chamberlength = None, fuelLength = None, preCombustionLength = None, portDiameter = None, nozzleType = None, inletRadius = None, throatRadius = None, areaRatio = None, outletAngle = None, exhaustk = None, exhaustMolecularWeight = None, efficiency = 0.9, thrustCurve = None, verbose = False):
            #Check the Motor Type and Find its Key
            if motorType in self.motorTypes:
                for key in range(len(self.motorTypes)):
                    if self.motorTypes[key] == motorType:
                        self.motorType = motorType
                        self.key = key
                        break
            else:
                print('Motor Type not Supported')
                print('Supported Motor Types: ' + str(self.motorTypes))
                self.error = True
                return
            
            #Set a Default Nozzle Inlet Radius if not Provided
            if not inletRadius and helper.validateNumber(type(chamberID)):
                if verbose == True:
                    print('Inlet radius not provided: using chamberID/2')
                inletRadius = chamberID/2
            
            #Check for a Valid Thrust Curve
            if type(thrustCurve) == np.ndarray:
                self.error = False
                if np.shape(thrustCurve)[0] != 2 and len(np.shape(thrustCurve)) < 10:
                    print('thrustCurve must be a 2 Dimensional array with time steps in row 1 and thrust in row 2. Thurst data over time should at least be 10 entries.')
                    self.error = True
                    return
                self.thrustCurve = thrustCurve
                self.index = 0
            elif motorType == 'commercial solid':
                print('A thrust curve must be provided for commercial solid motors')
                self.error = True
                return
            else:
                #Assign All Motor Attributes
                self.thrustCurve = None
                self.error = False
                self.on = True
                if motorType == 'solid':
                    self.chamber = Rocket.Motor.Chamber(motorType, chamberID, chamberlength, fuelLength, portDiameter, None, verbose)
                    self.error = self.chamber.error
                elif motorType == 'hybrid':
                    self.oxTank = Rocket.Motor.Tank(oxidizer, oxVolume, oxMass, oxPressure, oxOutletDiameter, oxPressurant, LOxFraction, verbose)
                    self.chamber = Rocket.Motor.Chamber(motorType, chamberID, chamberlength, fuelLength, portDiameter, preCombustionLength, verbose)
                    self.error = (self.oxTank.error == True or self.chamber.error == True)
                elif motorType == 'liquid':
                    self.oxTank = Rocket.Motor.Tank(oxidizer, oxVolume, oxMass, oxPressure, oxOutletDiameter, oxPressurant, LOxFraction, verbose)
                    self.fuelTank = Rocket.Motor.Tank(fuel, fuelVolume, fuelMass, fuelPressure, fuelOutletDiameter, fuelPressurant, None, verbose)
                    self.chamber = Rocket.Motor.Chamber(motorType, chamberID, chamberlength, None, None, None, verbose)
                    self.error = (self.oxTank.error == True or self.fuelTank.error == True or self.chamber.error == True)
                else:
                    self.error = True
                """We should calculate k and M for the exhaust here rather than making the user input it"""
                self.nozzle = Rocket.Motor.Nozzle(nozzleType, inletRadius, throatRadius, areaRatio, outletAngle, exhaustk, exhaustMolecularWeight, efficiency, verbose)
                self.error = (self.error or self.nozzle.error)
        
        #Update the Motor
        def updateMotor(self, a, theta, Patm, t, h):
            #Use Linear Interpolation on a Given Thrust Curve
            if np.any(self.thrustCurve):
                while self.thrustCurve[0][self.index] <= t and self.index < len(self.thrustCurve[0]) - 1:
                    self.index += 1
                w = (t - self.thrustCurve[0][self.index - 1])/(self.thrustCurve[0][self.index] - self.thrustCurve[0][self.index - 1])
                F = w*self.thrustCurve[1][self.index - 1] + (1 - w)*self.thrustCurve[1][self.index]
                mDot = None
                return F, mDot
            
            #Solid Motor Model
            elif self.motorType == 'solid':
                P0, T0 = self.chamber.updateSolid(h)
                return self.nozzle.updateNozzle(P0, T0, Patm)
            
            #Hybrid Motor Model
            elif self.motorType == 'hybrid':
                ODot = self.oxTank.updateTank(a, theta)
                P0, T0 = self.chamber.updateHybrid(ODot, h)
                return self.nozzle.updateNozzle(P0, T0, Patm)
            
            #Liquid Motor Model
            elif self.motorType == 'liquid':
                ODot = self.oxTank.updateTank(a, theta)
                FDot = self.fuelTank.updateTank(a, theta)
                P0, T0 = self.chamber.updateLiquid(ODot, FDot, h)
                return self.nozzle.updateNozzle(P0, T0, Patm)
            
            #Error State
            else:
                self.error = True
                return 0, 0
        
    class Airframe:
        #Create an airframe Object
        def __init__(self, emptyMass, cg, cp, diameter, cd):
            pass
        
        #Update the Forces on the Rocket
        def updateAirframe(self):
            pass
    
    #Create a Rocket Object
    def __init__(self, name, h = 0.01, rkOrder = 4, save = True, load = False):
        #Check the Name
        if type(name) != str:
            print('Rocket Name Must be a String')
            return
        
        #Load Rocket Information (read input config files)
        if load == True:
            pass
        
        #Create the Rocket
        self.airframe = Rocket.Airframe()
        self.motor = Rocket.Motor()
        
        #Attributes for the Rocket's Position. Array indexing corresponds with derivatives
        self.x = np.zeros(rkOrder + 2)
        self.y = np.zeros(rkOrder + 2)
        self.z = np.zeros(rkOrder + 2)
        self.theta = np.zeros(rkOrder)
        self.position = np.zeros(4)
        self.velocity = np.zeros(4)
        self.acceleration = np.zeros(4)
        self.jerk = np.zeros(4)
        
        #Speed Logging Variables
        self.railSpeed = 0
        self.maxSpeed = 0
        self.apogeeSpeed = 0
        self.drogueSpeed = 0
        self.mainSpeed = 0
        
        #Attributes for the Rocket's Status
        self.distance = 0
        self.condition = 0
        self.apogee = False
        self.landed = False
        
        #Other Simulation Control Variables
        self.rk = max(rkOrder, 2)
        self.h = h
        
        #Save the Rocket
        self.save()
    
    #Save a Rocket
    def save(self):
        pass
    
    #Update the Rocket's Position/Orientation Using its Old Position and Updated Derivatives
    def updateRocket(self):
        for i in range(self.rk - 1):
            for j in range(i + 1, self.rk):
                self.x[i] += (self.x[j]*self.h**(j - i))/factorial(j - i)
                self.y[i] += (self.y[j]*self.h**(j - i))/factorial(j - i)
                self.z[i] += (self.z[j]*self.h**(j - i))/factorial(j - i)
                self.theta[i] += (self.theta[j]*self.h**(j - i))/factorial(j - i)
        self.position = np.vstack((self.position, np.array([self.x[0], self.y[0], self.z[0], self.theta[0]])))
        self.velocity = np.vstack((self.velocity, np.array([self.x[1], self.y[1], self.z[1], self.theta[1]])))
        self.acceleration = np.vstack((self.acceleration, np.array([self.x[2], self.y[2], self.z[2], self.theta[2]])))
        self.jerk = np.vstack((self.jerk, np.array([self.x[3], self.y[3], self.z[3], self.theta[3]])))
    
    #Simulate a Single Flight
    def simulate(self, windSpeed, windDirection, gustSpeed, temperature, railDirection, railAngle = 6, railLength = 5.45, altitude = 1300, latitude = 32.94205, longitude = -106.91548, h = 0.1, static = False, plot = True):
        #Check the Static Fire Control
        if type(static) not in {bool, str}:
            print("Static Fire command must be: 'True'/'vertical', 'horizontal', or 'False'")
            return
        
        #Check Other Inputs
        
        #Check if Simulation is for a Static Fire or a Launch
        i = 0
        if static in {True, 'vertical'}:
            #Vertical Static Fire
            while self.Motor.on == True:
                self.Motor.updateMotor()
                i += 1
            
            #Show Stats Here
        elif static == 'horizontal':
            #Horizontal Static Fire
            pass
        else:
            #Simulate a Complete FLight
            
            #Set the Starting Position/Angle
            """I think we should make the launch rail the origin and then work with the GPS coordinates later
            self.x[0] = longitude
            self.y[0] = latitude
            self.z[0] = altitude
            """
            self.theta[0] = railAngle
            self.position = np.array([longitude, latitude, altitude, railAngle])
            
            #Simulate the Launch off the Rail
            while self.distance <= railLength:
                self.motor.updateMotor()
                self.updateAtmosphere()
                self.airframe.updateAirframe()
                self.updateRocket()
                i += 1
            
            #Simulate the Ascent Power
            while self.apogee == False:
                self.motor.updateRocket()
                self.updateAtmosphere()
                self.airframe.updateAirframe()
                self.updateRocket()
                i += 1
            
            #Simulate the Descent
            while self.landed == False:
                self.updateAtmosphere()
                self.airframe.updateAirframe()
                self.updateRocket
                i += 1
    
    #Simulate Several Flights
    def runSimulation(self, flights = 100, mapSize = (1000, 1000), plot = True):
        #Check the Number of Flights
        if type(flights) not in {int, np.int32} or flights < 1:
            print('The number of flights must be an integer greater than 0')
            return
        
        #Check the mapSize
        if type(mapSize) != tuple or len(mapSize) != 2 or type(mapSize[0]) not in {int, np.int32} or type(mapSize[0]) not in {int, np.int32} or mapSize[0] < 10 or mapSize[1] < 10:
            print('Map Size must be a 2-tuple with integer entries greater than 10')
            return
        
        #Initialize Array for Storing Sim Data
        #Rows are: Apogee, Rail Speed, Max Speed, Speed at Apogee, Speed Under Drogue, Speed Under Main, Landing Latitude, Landing Longitude
        results = np.zeros((flights, 8))
        
        
        #Simulate Several Flights
        for i in range(flights):
            results[i] = self.simulate(plot = False)
        
        #Calculate the Average Results
        averageResults = np.mean(results, axis = 0)
        
        #Create a Heat Map for the Landing Zone
        '''How can we position this apropriately and calculate the distance'''
        heatMap = np.zeros(mapSize)
        for i in range(mapSize[0]):
            for j in range(mapSize[1]):
                for k in range(flights):
                    #heatMap[i][j] += distance
                    pass
        
        #Normalize the Heat Map
        heatMap = max(heatMap) - heatMap
        heatMap /= np.sum(heatMap)
        
        #Color and Overlay the Heat Map