import numpy as np


def validateNumber(val):
    # add more types
    if type(val) in {int, float, np.int32, np.float64} and val >= 0:
        return True
    else:
        return False