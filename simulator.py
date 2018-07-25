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
from simulation import Simulation

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

            #Store variables
            self.simulation_runs = int(json_data["simulation_runs"])
            self.simulation_type = self.validateSimType(json_data["simulation_type"])
            self.simulation_output_dir = json_data["simulation_output_dir"]
            self.simulation_models = json_data["simulation_models"]
            self.rocket_specs = json_data["rocket_specs"]

            # STEP 2: run (instantiate) simulations according to the configs
            Simulation(self.simulation_runs,self.simulation_type,self.simulation_output_dir,self.simulation_models,self.rocket_specs)

        except:
            print(sys.exc_info())


    # HELPER FUNCTIONS
    def validateSimType(self, type):
        if type == "full_flight" or type == "hot_fire_test":
            return type
        else:
            print("invalid simulation type. Please choose one from 'fullFlight' and 'hotFireTest'")
            return None


def main():
    # sim = Simulator("simulator.json")
    #implement read file name from cmd line
    if len(sys.argv) == 2 and str(sys.argv[1]).endswith(".json"):
        print('Running simulator with config file:', sys.argv[1])
        sim = Simulator(sys.argv[1])
    else:
        print('Running simulator with default config file: simulator.json')
        sim = Simulator("simulator.json")


if __name__ == "__main__":
    main()
