# -*- coding: utf-8 -*-
"""
McGill Rocket Team Flight Simulator
Authors: Matt Saathoff, Mei Qi Tang

The kinematics model class describes a specific rocket's behaviour in specific environments. This class will contain
the GPS position, acceleration, rotation of a specific rocket in time.
"""

import numpy as np
# import matplotlib.pyplot as plt
# import IPython.display as IPy
# from math import factorial
# import helper
import json
from pprint import pprint
import sys

class KinematicsModel:

    # RocketState object
    # Position, rotation, momentum, and angular momentum
    # are passed to the rocket state class of the kinematics model from the rocket object
    class RocketState:
        def __init__(self, rocket):
            # TODO: validate inputs
            # Inputs should be numpy arrays
            self.X = rocket.position
            self.Q = rocket.rotation
            self.P = rocket.linearMomentum
            self.L = rocket.angularMomentum

    # RocketStateDerivative object
    # Linear velocity, angular velocity, force, and torque
    # KinematicsModel produces a RocketStateDerivative object
    class RocketStateDerivative:
        def __init__(self, rocketXDot, rocketQDot, rocketForce, rocketTorque):
            # TODO: Validate types
            # Outputs should be numpy arrays
            self.rocketXDot = rocketXDot
            self.rocketQDot = rocketQDot
            self.rocketForce = rocketForce
            self.rocketTorque = rocketTorque

    #STEP 1: read a rocket object, read time (normally t = 0)
    # Instantiate a KinematicsModel object
    def __init__(self, rocket, timeStep, inputFile="kinematics.json"):
        # TODO: Include parameters from kinematics.json
        try:
            with open("template_configs/" + inputFile) as json_file:
                json_data = json.load(json_file)

            pprint(json_data)

            self.rocket = rocket
            self.timeStep = timeStep
            self.model_constants = json_data["model_constants"]
            self.model_parameters = json_data["model_parameters"]
            self.rocketState = KinematicsModel.RocketState(self.rocket)
            self.RocketStateDerivative = KinematicsModel.model()

        except IOError:
            print("Error reading input config files when instantiating KinematicsModel object", sys.exc_info())

    # Helper methods
    def getTimeStepIndex(self):
        # get index value associated with current timeStep
        pass

    #STEP 2: compute kinematics
    # Function that computes the kinematic model and return a RocketStateDerivative object
    def model(self):
        # TODO insert logic for computing dynamics

        # timeStep index
        timeStepIndex = self.getTimeStepIndex()
        # State vectors
        X  = self.rocketState.X # [x, y, z]
        Q = self.rocketState.Q  # [s, v_x, v_y, v_z]
        s = Q[0]                # [s]
        v = Q[1:]               # [v_x, v_y, v_z]
        P = self.rocketState.P  # [P_x, P_y, P_z]
        L = self.rocketState.L  # [L_x, L_y, L_z]

        # Refernce Yaw, Pitch, Roll axis vectors
        yaw_0 = self.model_constants["referenceYawAxis"]
        pitch_0 = self.model_constants["referencePitchAxis"]
        roll_0 = self.model_constants["referenceRollAxis"]

        # Model parameters that varry with time
        # TODO: Need to access value corresponding to current timeStep
        mass = self.model_parameters["mass"][timeStepIndex]
        yaw_Ixx = self.model_parameters["yawAxisMomentOfInertia"][timeStepIndex]
        pitch_Iyy = self.model_parameters["pitchAxisMomentOfInertia"][timeStepIndex]
        roll_Izz = self.model_parameters["rollAxisMomentOfInertia"][timeStepIndex]

        # Rotation matrix
        rotationMatrix = np.zeros((3,3))
        rotationMatrix[0,0] = 1 - 2*(Q[2]**2 - Q[3]**2) # row 1, col 1
        rotationMatrix[0,1] = 2*(Q[1]*Q[2] - Q[0]*Q[3]) # row 1, col 2
        rotationMatrix[0,2] = 2*(Q[1]*Q[3] + Q[0]*Q[2]) # row 1, col 3
        rotationMatrix[1,0] = 2*(Q[1]*Q[2] + Q[0]*Q[3]) # row 2, col 1
        rotationMatrix[1,1] = 1 - 2*(Q[1]**2 - Q[3]**2) # row 2, col 2
        rotationMatrix[1,2] = 2*(Q[2]*Q[3] - Q[0]*Q[1]) # row 2, col 3
        rotationMatrix[2,0] = 2*(Q[1]*Q[3] - Q[0]*Q[2]) # row 3, col 1
        rotationMatrix[2,1] = 2*(Q[2]*Q[3] + Q[0]*Q[1]) # row 3, col 2
        rotationMatrix[2,2] = 1 - 2*(Q[1]**2 - Q[2]**2) # row 2, col 3

        # Reference inertia tensor
        refIntertiaTensor = np.zeros((3,3))
        refIntertiaTensor[0,0] = yaw_Ixx
        refIntertiaTensor[1,1] = pitch_Iyy
        refIntertiaTensor[2,2] = roll_Izz

        # Updated Yaw, Pitch, Roll axis vectors
        yaw = np.matmul(rotationMatrix, np.transpose(yaw_0))
        pitch = np.matmul(rotationMatrix, np.transpose(pitch_0))
        roll = np.matmul(rotationMatrix, np.transpose(roll_0))

        # Linear velocity
        xDot = P / mass

        # Angular velocity
        omega = np.matmul(np.matmul(rotationMatrix, np.linalg.inv(refIntertiaTensor)), np.transpose(np.matmul(L, rotationMatrix)))
        sDot = 0.5*np.dot(omega, v)
        vDot = 0.5*(s*omega + np.cross(omega, v))
        qDot = np.append(sDot, vDot)

        # Angle of attack

        # Reynolds number

        # Mach number

        # Fin angle of attack

        # Force
        force = 0

        # Torque
        torque = 0

    #STEP 3: return GPS, accel, rot of a rocket at any instant in time using the environment of that moment
        return KinematicsModel.RocketStateDerivative(xDot, qDot, force, torque)
