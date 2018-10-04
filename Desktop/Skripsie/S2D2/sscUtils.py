import numpy as np
import mappingUtils
import matplotlib.pyplot as plt

def getOptimalConeRad(forceConeRad):
    # Assume that the cone rad falls between -PI and PI
    tanFCR = np.tan(forceConeRad)
    if (np.abs(tanFCR) < 1e-6):
        OCR = 0.0
    else:
        tanFCR2 = tanFCR**2.0
        topTerm = -3.0 + np.sqrt(9.0 + 8.0*tanFCR2)
        botTerm = 4 * tanFCR
        tanOCR = topTerm / botTerm
        OCR = np.arctan(tanOCR)
    return OCR

def coneRadForMaxF2D(rayDir, fDir):
    cosOfForceConeRad = fDir.dot(rayDir) #By reducing sunDir magnitude to 1 in previous step, no division is necessary
    if (cosOfForceConeRad > 1.0):
        forceConeRad = 0.0
    elif (cosOfForceConeRad < -1.0):
        forceConeRad = np.pi
    else:
        forceConeRad = np.arccos(cosOfForceConeRad)
    if (forceConeRad > np.pi/2): #Can't be obtained with solar sail
        forceConeRad = (np.pi/2 - 1e-6)
    elif (forceConeRad < (-np.pi/2)):
        forceConeRad = (-np.pi/2 + 1e-6)
    #print(forceConeRad)
    OCR = getOptimalConeRad(forceConeRad)
    return OCR

def coneRadToDir2D(rayDir, coneRad):
    rayRad = mappingUtils.dirToRadians2D(rayDir)
    coneDir = mappingUtils.radiansToDir2D(rayRad + coneRad)
    return coneDir

def testOCR(numSteps):
    startRad = -np.pi/2
    radChange = np.pi
    radChange = radChange / numSteps
    radArray = np.zeros(numSteps)
    ocrArray = np.zeros(numSteps)
    for i in range(numSteps):
        rad = startRad + i*radChange
        radArray[i] = rad * 180 / np.pi
        ocrArray[i] = getOptimalConeRad(rad) * (180/np.pi)

    plt.plot(radArray, ocrArray)
    plt.show()

def testEccOCR(numSteps):
    startRad = -2 * np.pi
    radChange = 1.5 * np.pi
    radChange = radChange / numSteps
    radArray = np.zeros(numSteps)

    rDirArray = np.zeros((numSteps, 2))
    fDirArray = np.zeros((numSteps, 2))

    rAngleArray = np.zeros(numSteps)
    fAngleArray = np.zeros(numSteps)
    dotArray = np.zeros(numSteps)
    coneArray = np.zeros(numSteps)

    for i in range(numSteps):
        rad = startRad + i * radChange
        rDir = mappingUtils.radiansToDir2D(rad)
        fDir = mappingUtils.radiansToDir2D(rad + np.pi/2)
        rAngle = mappingUtils.dirToRadians2D(rDir) * (180/np.pi)
        fAngle = mappingUtils.dirToRadians2D(fDir) * (180/np.pi)
        dP = rDir.dot(fDir)
        cA = coneRadForMaxF2D(rDir, fDir)

        radArray[i] = rad
        rDirArray[i] = rDir
        fDirArray[i] = fDir
        rAngleArray[i] = rAngle
        fAngleArray[i] = fAngle
        dotArray[i] = dP * 100
        coneArray[i] = cA * (180/np.pi)

    pltR = rDirArray.T
    pltF = fDirArray.T * 1.2
    '''
    plt.plot(rAngleArray)
    plt.plot(fAngleArray)
    plt.plot(dotArray)
    plt.show()
    '''
    plt.plot(coneArray)
    plt.show()

#testOCR(500)
#def testOptAngleSmalls(numSteps):
#testEccOCR(360)
