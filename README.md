# simulator
### Code base for McGill Rocket Team's flight simulator

#### basic code structure

1. run program on cmdline "python simulator.py"
2. simulator.py instantiates simulations from simulation.py
3. simulation.py instantiates Rocket object and models (eg. kinematic and atmospheric models)

```
Flow: Simulator() <--(reads)  simulator.json
            |      \
            |       \        ... ( x simulation_runs)
            v        v
      Simulation()  Simulation()   ...
       /    |
      /     |
     v      v
Rocket()  Models: - KinematicModel() <--(reads)  kinematics.json
                  - AtmosphericModel() <--(reads)  atmospherics.json
                  - ...

```
Each of these classes reads JSON format config files. They look for default ones in the template_configs folder.
If you want to use different config files, you can:
1. run the program with the simulator config file as a cmdline argument like "python simulator.py myConfig.json"
2. in the simulator config file that you decide to use, specify the other config files that you want. Eg. "simulation_models": {"kinematic":"myKinematics.json", "atmospheric":"myAtmospherics.json"}

