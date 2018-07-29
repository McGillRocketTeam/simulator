import numpy as np

class Helper:
    def __init__(self):
        self.motorTypes = ['commercial solid', 'solid', 'hybrid', 'liquid']
        self.liquidFuels = ['RP1', 'ethanol', 'methanol', 'LNG', 'gasoline', 'N2H4', 'MMH', 'UDMH']
        self.hybridFuels = ['paraffin', 'HTPB', 'PE']
        self.oxidizers = ['N2O', 'LOx', 'nitrox', 'H2O2', 'ClF3', 'ClF5', 'N2O4', 'N2F4', 'ClO3F', 'FClO4', 'FLOx']
        self.pressurants = ['He', 'N2', 'Ar']
        self.nozzleTypes = ['CD']
        self.R = 8.3144598  # J/K*mol

    def containsValue(self, list, element):
        try:
            if list == 'motorTypes' and element in self.motorTypes:
                return True
            elif list == 'liquidFuels'and element in self.liquidFuels:
                return True
            elif list == 'hybridFuels'and element in self.hybridFuels:
                return True
            elif list == 'oxidizers'and element in self.oxidizers:
                return True
            elif list == 'pressurants'and element in self.pressurants:
                return True
            elif list == 'nozzleTypes'and element in self.nozzleTypes:
                return True
            else:
                return False
        except:
            print("Error encountered while checking is element is in list")
            return False

    def getIndex(self, list, element):
        try:
            if list == 'motorTypes':
                return self.motorTypes.index(element)
            elif list == 'liquidFuels':
                return self.liquidFuels.index(element)
            elif list == 'hybridFuels':
                return self.hybridFuels.index(element)
            elif list == 'oxidizers':
                return self.oxidizers.index(element)
            elif list == 'pressurants':
                return self.pressurants.index(element)
            elif list == 'nozzleTypes':
                return self.nozzleTypes.index(element)
            else:
                return False
        except ValueError:
            print("Can't get index of an element not in list")
            return False
        except:
            print("Error encountered while getting index of element in list")
            return False

    def getConstant(self, constant):
        if constant == 'R':
            return self.R
        else:
            print('Constant "' + constant + '" does not exist')
            return False

    def validateNumber(self, val):
        # add more types
        if type(val) in {int, float, np.int32, np.float64} and val >= 0:
            return True
        else:
            return False