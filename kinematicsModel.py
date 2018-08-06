# -*- coding: utf-8 -*-
"""
McGill Rocket Team Flight Simulator
Authors: Matt Saathoff, Mei Qi Tang, Luca D'Angelo

The kinematics model class describes a specific rocket's behaviour in specific environments. This class will contain
the GPS position, acceleration, rotation of a specific rocket in time.
"""

import numpy as np
# import matplotlib.pyplot as plt
# import IPython.display as IPy
# from math import factorial
# import helper

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
            self.rocketXDot = rocketVelocity
            self.rocketQDot = rocketQuaternionDerivative
            self.rocketForce = rocketForce
            self.rocketTorque = rocketTorque

#STEP 1: read a rocket object, read time (normally t = 0)
    # Instantiate a KinematicsModel object
    def __init__(self, rocket, time, inputFile="kinematics.json"):
        # TODO: Include parameters from kinematics.json
        self.rocket = rocket
        self.timeStep = time
        self.rocketState = KinematicsModel.RocketState(self.rocket)
        self.rocketStateDerivative = KinematicsModel.model()

    # Function that computes the kinematic model and return a RocketStateDerivative object
    def model(self):
        # TODO insert logic for computing dynamics
        return KinematicsModel.RocketStateDerivative(xDot, qDot, force, torque)

#STEP 2: compute kinematics

#STEP 3: return GPS, accel, rot of a rocket at any instant in time using the environment of that moment
