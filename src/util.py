import numpy as np
import numbers
import math

def make_string(val):
    """converts the value into a string if it is not already a string"""
    if val and not isinstance(val,basestring) :
        result = str(val)
    else:
        result = val

    return result


def isnan(val):
    """Multiple tests if it is nan, either math or nump"""
    return  isinstance(val,numbers.Number) and (np.isnan(val) or math.isnan(val))
