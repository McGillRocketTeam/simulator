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
import os
import datetime
from rocket import Rocket
from kinematicsModel import KinematicsModel
from atmosphericModel import AtmosphericModel


class Simulation:

    # STEP 1: create output dir & sub folder to store simulation results
    def __init__(self, runs, type, outputDir, models, rocket):
        dir_name= "output_results/" + datetime.datetime.now().strftime("%Y-%m-%d_%H.%M.%S") #using '.' since dir name doesn't allow to contain ':'
        if outputDir:
            dir_name = "output_results/" + outputDir

        os.mkdir(dir_name)

        #STEP 2: simulate runs times and store flight results to outputdir/simulation_{#}.csv
        for i in range(0, runs):
            f = open(dir_name + '/simulation_' + str(i) + '.csv', 'w')         #async?
            f.write(self.simulate(type, models, rocket))


    #STEP 3: instantiate according to configs and link desired models together
    def simulate(self, type, models, rocket):
        if rocket:
            self.rocket = Rocket(rocket)
        if models["kinematic"]:
            self.kinematicsModel = KinematicsModel("")
        if models["atmospheric"]:
            self.atmosphericModel = AtmosphericModel("")

        #TODO: implement simulator logic here

        return ""  #return simulation variables here
