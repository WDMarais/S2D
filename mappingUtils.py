from scipy.constants import au
from miscConstants import DIMENSION_ERROR
from miscConstants import LAZY_CODER_ERROR

import numpy as np

import sys
from miscConstants import INVALID_PARAMETER_ERROR

def scaleDistFromSI(newUnit):
    if (newUnit == 'm'):
        scaleFactor = 1.0
    elif (newUnit == 'au'):
        scaleFactor = (1.0 / au)
    elif (newUnit == 'km'):
        scaleFactor = (1e-3)
    else:
        print("mappingUtils scaleDistFromSI error - invalid parameter")
        print(newUnit)
        sys.exit(INVALID_PARAMETER_ERROR)
    return scaleFactor

def scaleTimeFromSI(newUnit):
    if (newUnit == 'sec'):
        scaleFactor = 1.0
    elif (newUnit == 'min'):
        scaleFactor = (1.0/60.0)
    elif (newUnit == 'hr'):
        scaleFactor = (1.0/(60.0 * 60.0))
    elif (newUnit == 'day'):
        scaleFactor = (1.0/(60.0 * 60.0 * 24.0))
    elif (newUnit == 'mon'):
        scaleFactor = (1.0/(60.0 * 60.0 * 24.0 * 30))
    elif (newUnit == 'yr'):
        scaleFactor = (1.0/(60.0 * 60.0 * 24.0 * 365.25))
    elif (newUnit == 'cty'):
        scaleFactor = (1.0/(60.0 * 60.0 * 24.0 * 365.25 * 100))
    else:
        print("mappingUtils scaleTimeFromSI error - invalid parameter")
        print(newUnit)
        sys.exit(INVALID_PARAMETER_ERROR)
    return scaleFactor

def scaleVelFromSI(newUnit):
    distUnit, timeUnit = newUnit.split('/')
    distScale = scaleDistFromSI(distUnit) # Dist Units per meter
    timeScale = scaleTimeFromSI(timeUnit) # Time Units per second
    velScaleFactor = distScale / timeScale
    return velScaleFactor

def scaleAngleFromRad(newUnit):
    scaleFactor = 1.0
    if (newUnit == 'deg'):
        scaleFactor = (180/0 / np.pi)
    elif (newUnit == 'rev'):
        scaleFactor = (1.0/(2.0*np.pi))
    return scaleFactor

def modRad(radians):
    return (radians % (2*np.pi))

def radModPlusMinusPi(radians):
    modded = modRad(radians)
    if (modded > np.pi):
        modded = modded - 2 * np.pi
    return modded

def dirToAngle(rArray):
    dim = len(rArray)
    if (dim == 2):
        x = rArray[0]
        y = rArray[1]
        angle = np.arctan2(y, x)
        return angle
    else:
        dimStr = str(dim) + "D"
        print("------------------------------")
        print("mappingUtils dirToAngle Error")
        print("Was expecting 2D rArray, received", dimstr, "instead")
        print("-------------------------------")
        sys.exit(LAZY_CODER_ERROR)

def dirToRadians2D(rArray):
    x = rArray[0]
    y = rArray[1]
    radians = np.arctan2(y, x)
    return radians

def radiansToDir2D(radians):
    x = np.cos(radians)
    y = np.sin(radians)
    rDir = np.array([x, y])
    return rDir

def rotate2D(vector, radians):
    x, y = vector
    c, s = np.cos(radians), np.sin(radians)
    rx = x * c - y * s
    ry = x * s + y * c
    return np.array([rx, ry])
