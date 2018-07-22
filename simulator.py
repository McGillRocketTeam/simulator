# -*- coding: utf-8 -*-
"""
McGill Rocket Team Flight Simulator
Authors: Matt Saathoff, Mei Qi Tang

The simulator class configures the simulations. After reading input configuration instructions, it will call the
simulation class to execute simulations according to the configs. Possible configurable variables could be the number
of simulations to run (eg. 100), the specific modules active in each simulation (eg. atmospheric model), the type of simulation
to run (eg. vertical fire test), etc.
"""

# import numpy as np
# import matplotlib.pyplot as plt
# import IPython.display as IPy
# from math import factorial
# import helper
import json
from pprint import pprint
import sys

class Simulator:
    def __init__(self, inputFile):
    #STEP 1: read, validate, and store values from input JSON config file
        try:
            if inputFile:
                with open(inputFile) as json_file:
                    json_data = json.load(json_file)
            else:
                with open("simulator.json") as json_file:
                    json_data = json.load(json_file)

            pprint(json_data)                       #pretty-print
            print(json_data["rocket_specs"])

            #Store variables
            self.simulation_runs = int(json_data["simulation_runs"])
            self.simulation_type = self.validateSimType(json_data["simulation_type"])
            self.simulation_output_dir = json_data["simulation_output_dir"]
            self.simulation_modules = json_data["simulation_modules"]
            self.rocket_specs = json_data["rocket_specs"]

        except:
            print(sys.exc_info())

    #STEP 2: run (instantiate) simulations according to the configs

    #Helper Functions

    def validateSimType(self, type):
        if type == "fullFlight" or type == "hotFireTest":
            return type
        else:
            print("invalid simulation type. Please choose one from 'fullFlight' and 'hotFireTest'")
            return None


def main():
    sim = Simulator("simulator.json")

if __name__ == "__main__":
    main()
