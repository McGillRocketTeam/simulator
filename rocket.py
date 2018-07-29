# -*- coding: utf-8 -*-
"""
McGill Rocket Team Flight Simulator
Authors: Matt Saathoff, Mei Qi Tang

The Rocket class creates a rocket model that includes a Motor and Airframe. The Motor is further separated
in a Tank, Chamber, and Nozzle. This class only contains the hard "specs" about the rocket. It only defines
the properties of the rocket and not its flight behaviour for example.
"""

import numpy as np
# import matplotlib.pyplot as plt
# import IPython.display as IPy
# from math import factorial
from helper import Helper
from math import factorial
import json
from pprint import pprint
import sys

helper = Helper()

# A Class for Creating Rocket Objects
class Rocket:
    # A Class for Creating Motor Objects
    class Motor:
        # A Class for Creating Tank Objects
        class Tank:
            # Create a tank Object
            def __init__(self, motorType, configObject, verbose=False):
                # Error State Attribute
                self.error = False
                pressurant = None
                LOxFraction = None
                propellant = configObject["oxName"]
                volume = configObject["oxVolume"]
                propellantMass = configObject["oxMass"]
                initialPressure = configObject["oxPressure"]
                outletDiameter = configObject["oxOutletDiameter"]
                if configObject["oxPressurant"]:
                    pressurant = configObject["oxPressurant"]
                if configObject["LOxFraction"]:
                    LOxFraction = configObject["LOxFraction"]

                # Check/Assign the Propellant Type/Key
                if motorType == 'liquid' and helper.containsValue('liquidFuels', propellant):
                    self.type = 'liquid_fuel'
                    self.propellant = propellant
                    self.key = helper.getIndex('liquidFuels', propellant)
                elif motorType == 'hybrid' and helper.containsValue('hybridFuels', propellant):
                    self.type = 'hybrid_fuel'
                    self.propellant = propellant
                    self.key = helper.getIndex('hybridFuels', propellant)
                elif motorType in {'liquid', 'hybrid'}:
                    self.error = True
                    if motorType == 'liquid':
                        print('Invalid Fuel Selection\nSupported Liquid Fuels: ' + str(helper.liquidFuels))
                    else:
                        print('Invalid Fuel Selection\nSupported Hybrid Fuels: ' + str(helper.hybridFuels))
                    return
                elif helper.containsValue('oxidizers', propellant):
                    self.type = 'oxidizer'
                    self.propellant = propellant
                    self.key = helper.getIndex('oxidizers', propellant)
                    if propellant in {'nitrox', 'FLOx'}:
                        if type(LOxFraction) not in {float, np.float64} or LOxFraction < 0 or LOxFraction > 1:
                            print('LOxFraction must be a number between 0 and 1')
                            self.error = True
                            return
                else:
                    if motorType == 'liquid':
                        print('Propellants Not Supported\nSupported Liquid Fuels: ' + str(helper.liquidFuels))
                    else:
                        print('Propellants Not Supported\nSupported Hybrid Fuels: ' + str(helper.hybridFuels))
                    print('Supported Oxidizers: ' + str(helper.oxidizers))
                    self.error = True
                    return

                # Check the Volume
                if helper.validateNumber(volume):
                    self.V = volume
                else:
                    print('Volume must be a number greater than 0')
                    self.error = True
                    return

                # Check/Assign the Propellant Mass
                if helper.validateNumber(propellantMass):
                    self.m = np.array([propellantMass])
                else:
                    print('Propellent Mass must be a number greater than 0')
                    self.error = True
                    return

                # Check/Assign the Initial Pressure
                if helper.validateNumber(initialPressure):
                    self.P = np.array([initialPressure])
                else:
                    print('Initial Pressure must be a number greater than 0')
                    self.error = True
                    return

                # Check/Assign the Outlet Diameter
                if helper.validateNumber(outletDiameter):
                    self.D = outletDiameter
                else:
                    print('Outlet Diameter must be a number greater than 0')
                    self.error = True
                    return

                # Check/Assign the Pressurant
                if helper.containsValue('pressurants', pressurant):
                    self.pressurant = pressurant
                else:
                    if pressurant:
                        print('Pressurant Not Supported')
                        print('Supported Pressurants: ' + str(helper.pressurants))
                    else:
                        print("No inputted pressurant")
                    self.error = True

                # Assign the Verbose Attribute
                self.verbose = verbose

            # Initialize a Tank Before Launch
            def initialize(self, temperature):
                pass

            # Update a Tank
            def updateTank(self, a, theta):
                pass

        # A Class for Creating Combustion Chamber Objects
        class Chamber:
            # Create a chamber Object
            def __init__(self, motorType, configObject, verbose=False):
                # Error State Attribute
                self.error = False
                fuelLength = None
                portDiameter = None
                preCombustion = None
                ID = configObject["chamberID"]
                length = configObject["chamberLength"]
                self.verbose = verbose

                #TODO: validate variables and assign to self


            # Solid Motor Combustion Chamber Model
            def updateSolid(self, h):
                return None, None

            # Hybrid Motor Combustion Chamber Model
            def updateHybrid(self, ODot, h):
                return None, None

            # Liquid Motor Combustion Chamber Model
            def updateLiquid(self, ODot, FDot, h):
                return None, None

        # A Class for Creating Nozzle Objects
        class Nozzle:
            # Create a Nozzle Object
            def __init__(self, configObject, verbose=False):
                efficiency = 1
                nozzleType = configObject["nozzleType"]
                inletRadius = configObject["inletRadius"]
                throatRadius = configObject["throatRadius"]
                areaRatio = configObject["areaRatio"]
                outletAngle = configObject["outletAngle"]
                k = configObject["exhaustk"]
                molecularWeight = configObject["exhaustMolecularWeight"]
                if configObject["efficiency"]:
                   efficiency = configObject["efficiency"]

                # Check the Nozzle Type and get its Key
                if helper.containsValue('nozzleTypes', nozzleType):
                    self.key = helper.getIndex('nozzleTypes', nozzleType)
                else:
                    print('Nozzle Type Not Supported')
                    print('Supported Nozzle Types: ' + str(helper.nozzleTypes))
                    self.error = True
                    return

                # Check the Throat Radius
                if not helper.validateNumber(throatRadius):
                    print('Throat Radius must be a number greater than 0')
                    self.error = True
                    return

                # Check the Inlet Radius
                if not helper.validateNumber(inletRadius) or inletRadius <= throatRadius:
                    print('Inlet Radius must be a number greater than the Throat Radius')
                    self.error = True
                    return

                # Check the Area Ratio
                if not helper.validateNumber(areaRatio) or areaRatio <= 1:
                    print('Area Ratio must be a number greater than 1')
                    self.error = True
                    return

                # Check the Outlet Angle
                if not helper.validateNumber(outletAngle) or outletAngle < 0 or outletAngle > 90:
                    print('Outlet Angle must be a number between 0 and 90 degrees')
                    self.error = True
                    return

                # Check k
                if not helper.validateNumber(k) or k < 1:
                    print('k must be a number greater than 1')
                    self.error = True
                    return

                # Check the Molecular Weight
                if not helper.validateNumber(molecularWeight) or molecularWeight < 0:
                    print('Molecular weight must be a number greater than 0')
                    self.error = True
                    return

                # Check the Efficiancy
                if not helper.validateNumber(efficiency) or efficiency > 1:
                    print('Efficiency must be a number between 0 and 1')
                    self.error = True
                    return

                # Assign All Attributes
                self.nozzleType = nozzleType
                self.epsilon = areaRatio
                self.ri = inletRadius
                self.rt = throatRadius
                self.re = throatRadius * np.sqrt(areaRatio)
                self.theta = outletAngle
                self.k = k
                self.M = molecularWeight
                self.R = helper.getConstant("R")
                self.F = 0
                self.mDot = 0
                self.efficiency = efficiency
                self.choked = False
                self.error = False
                self.verbose = verbose

            # Update the Nozzle (Assuming Isentropic Flow)
            def updateNozzle(self, P0, T0, Patm, precision=10 ** -10):
                # Check for Choked Flow
                self.choked = (Patm <= P0 * (2 / (self.k + 1)) ** ((self.k + 1) / (self.k - 1)))

                # Determine the Thrust and Mass Flow
                if self.choked:
                    # Calculate the Mass FLow through the Throat
                    self.mDot = self.efficiency * P0 * np.pi * (self.rt ** 2) * np.sqrt(self.k / (self.R * T0)) * (
                                1 + (self.k - 1) / 2) ** ((self.k + 1) / (2 - 2 * self.k))

                    # Calculate the Exit Mach
                    """
                    There is almost certainly a better way to solve this, but this method isn't bad
                    either since it runs in O(log(Me/precision))
                    """
                    Max = 1
                    while self.epsilon > (
                            (1 + 0.5 * (self.k - 1) * Max ** 2) ** (0.5 * (self.k + 1) / (self.k - 1))) / Max:
                        Max *= 2
                    Min = Max / 2
                    while abs(Max - Min) > precision:
                        Me = (Min + Max) / 2
                        if self.epsilon > (
                                (1 + 0.5 * (self.k - 1) * Me ** 2) ** (0.5 * (self.k + 1) / (self.k - 1))) / Me:
                            Min = Me
                        else:
                            Max = Me

                    Me = (Min + Max) / 2  # TODO: not sure about this here...
                    # Calculate Pe/P0 (Exhaust Pressure over Stagnation Pressure)
                    Pr = ((2 / (self.k + 1)) / (1 + 0.5 * (self.k - 1) * Me ** 2)) ** (self.k / (self.k - 1))

                    # Calculate the Exhaust Velocity
                    Ve = np.sqrt(
                        2 * self.R * T0 * self.k * ((1 - (Pr) ** ((self.k - 1)) / self.k)) / (self.M * (self.k - 1)))

                    # Calculate the Thrust
                    self.F = self.mDot * Ve + (P0 * Pr - Patm) * np.pi * self.re ** 2
                else:
                    pass

        # Create a Motor Object
        def __init__(self, configObject):

            motorType = configObject["motorType"]
            thrustCurve = configObject["thrustCurve"]
            self.fuel_config = configObject["fuel"]
            self.oxidizer_config = configObject["oxidizer"]
            self.chamber_config = configObject["chamber"]
            self.nozzle_config = configObject["nozzle"]
            self.verbose = configObject["verbose"]

            # Check the Motor Type and Find its Key
            if helper.containsValue('motorTypes', motorType):
                self.motorType = motorType
                self.key = helper.getIndex('motorTypes', motorType)
            else:
                print('Motor Type not Supported')
                print('Supported Motor Types: ' + str(helper.motorTypes))
                self.error = True
                return

            # Check for a Valid Thrust Curve
            if type(thrustCurve) == np.ndarray:
                self.error = False
                if np.shape(thrustCurve)[0] != 2 and len(np.shape(thrustCurve)) < 10:
                    print(
                        'thrustCurve must be a 2 Dimensional array with time steps in row 1 and thrust in row 2. Thurst data over time should at least be 10 entries.')
                    self.error = True
                    return
                self.thrustCurve = thrustCurve
                self.index = 0
            elif motorType == 'commercial solid':
                print('A thrust curve must be provided for commercial solid motors')
                self.error = True
                return
            else:
                # Assign All Motor Attributes
                self.thrustCurve = None
                self.error = False
                self.on = True
                if motorType == 'solid':
                    self.chamber_config["preCombustionLength"] = None
                    self.chamber = Rocket.Motor.Chamber(motorType, self.chamber_config, self.verbose)
                    self.error = self.chamber.error
                elif motorType == 'hybrid':
                    self.oxTank = Rocket.Motor.Tank(self.oxidizer_config, self.verbose)
                    self.chamber = Rocket.Motor.Chamber(motorType, self.chamber_config, self.verbose)
                    self.error = (self.oxTank.error == True or self.chamber.error == True)
                elif motorType == 'liquid':
                    self.fuel_config["pressurant"] = None
                    self.chamber_config["fuelLength"] = None
                    self.chamber_config["preCombustionLength"] = None
                    self.chamber_config["portDiameter"] = None
                    self.oxTank = Rocket.Motor.Tank(self.oxidizer_config, self.verbose)
                    self.fuelTank = Rocket.Motor.Tank(self.fuel_config, self.verbose)
                    self.chamber = Rocket.Motor.Chamber(motorType, self.chamber_config, self.verbose)
                    self.error = (self.oxTank.error == True or self.fuelTank.error == True or self.chamber.error == True)
                else:
                    self.error = True

                # Set a Default Nozzle Inlet Radius if not Provided
                if not self.nozzle_config["inletRadius"] and helper.validateNumber(type(self.chamber_config["chamberID"])):
                    if self.verbose == True:
                        print('Inlet radius not provided: using chamberID/2')
                        self.nozzle_config["inletRadius"] = self.chamber_config["chamberID"] / 2
                """We should calculate k and M for the exhaust here rather than making the user input it"""
                self.nozzle = Rocket.Motor.Nozzle(self.nozzle_config, self.verbose)
                self.error = (self.error or self.nozzle.error)

        # Update the Motor
        def updateMotor(self, a, theta, Patm, t, h):
            # Use Linear Interpolation on a Given Thrust Curve
            if np.any(self.thrustCurve):
                while self.thrustCurve[0][self.index] <= t and self.index < len(self.thrustCurve[0]) - 1:
                    self.index += 1
                w = (t - self.thrustCurve[0][self.index - 1]) / (
                            self.thrustCurve[0][self.index] - self.thrustCurve[0][self.index - 1])
                F = w * self.thrustCurve[1][self.index - 1] + (1 - w) * self.thrustCurve[1][self.index]
                mDot = None
                return F, mDot

            # Solid Motor Model
            elif self.motorType == 'solid':
                P0, T0 = self.chamber.updateSolid(h)
                return self.nozzle.updateNozzle(P0, T0, Patm)

            # Hybrid Motor Model
            elif self.motorType == 'hybrid':
                ODot = self.oxTank.updateTank(a, theta)
                P0, T0 = self.chamber.updateHybrid(ODot, h)
                return self.nozzle.updateNozzle(P0, T0, Patm)

            # Liquid Motor Model
            elif self.motorType == 'liquid':
                ODot = self.oxTank.updateTank(a, theta)
                FDot = self.fuelTank.updateTank(a, theta)
                P0, T0 = self.chamber.updateLiquid(ODot, FDot, h)
                return self.nozzle.updateNozzle(P0, T0, Patm)

            # Error State
            else:
                self.error = True
                return 0, 0

    class Airframe:
        # Create an airframe Object
        def __init__(self, configObject):
            # TODO: validate types
            self.emptyMass = float(configObject["emptyMass"])
            self.cg = float(configObject["cg"])
            self.cp = float(configObject["cp"])
            self.diameter = float(configObject["diameter"])
            self.cd = float(configObject["cd"])
            pass

        # Update the Forces on the Rocket
        def updateAirframe(self):
            pass

    # STEP 1: instantiate a Rocket Object according to input configs
    def __init__(self, inputFile="rocket.json"):
        try:
            with open("template_configs/" + inputFile) as json_file:
                json_data = json.load(json_file)

            pprint(json_data)

            self.motor_specs = json_data["motor_specs"]
            self.airframe_specs = json_data["airframe_specs"]
            self.name = str(json_data["name"])
            self.verbose = json_data["verbose"]

            # Check the Name
            if not self.name:
                print('Using default rocket name: "my_rocket"')
                self.name = 'my_rocket'

            self.airframe = Rocket.Airframe(self.airframe_specs)
            self.motor = Rocket.Motor(self.motor_specs)

        except IOError:
            print("Error reading input config files when instantiating Rocket object", sys.exc_info())

    # Save a Rocket
    def save(self):
        pass

    # Update the Rocket's Position/Orientation Using its Old Position and Updated Derivatives
    def updateRocket(self):
        pass

