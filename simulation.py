# -*- coding: utf-8 -*-
"""
McGill Rocket Team Flight Simulator
Authors: Matt Saathoff, Mei Qi Tang

Every instance of the simulation class is responsible to run a single simulation. The simulation class puts all
the models (eg. atmospheric and kinematic) together to generate snapshots of data in time using linear interpolation.
After simulation, it will return a final snapshot of the important flight variables (eg. apogee, gps coordinates,
flight time, both drogue and main deployment times, coordinates, etc).
"""

# import numpy as np
# import matplotlib.pyplot as plt
# import IPython.display as IPy
# from math import factorial
# import helper

#STEP 1: instantiate according to configs = instantiate and link the right models together

#STEP 2: return the variables we care about
