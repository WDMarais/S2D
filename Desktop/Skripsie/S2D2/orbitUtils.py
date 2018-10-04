import numpy as np
import datetime
from scipy.constants import au
import mappingUtils

def keplerElliptical(M, ecc, tol, maxIter):
    dE = 1
    E = M
    loops = 0
    while((np.abs(dE) > tol) and (loops < maxIter)):
        loops += 1
        dE = (E - ecc * np.sin(E) - M)/(1 - ecc * np.cos(E))
        E -= dE

    return E

def solveKepler(M, ecc, tol=1e-6, maxIter=5):
    if (ecc < 0.9):
        E = keplerElliptical(M, ecc, tol, maxIter)
    else:
        raise Exception("Eccentricity near-hyperbolic")
    return E

def addCorrectionFactor(oE, deltaCenturies, b, c, s, f):
    newOE = oE
    newOE += b * (deltaCenturies**2.0)
    newOE += c * np.cos(f * deltaCenturies)
    newOE += s * np.sin(f * deltaCenturies)
    return newOE

def calculateJPLOE(oE0, oEDot, deltaCenturies, cF=False):
    oE = oE0 + oEDot * deltaCenturies
    if not (cF == False):
        b, c, s, f = cF["b"], cF["c"], cF["s"], cF["f"]
        oE = addCorrectionFactor(b, c, s, f)
    return oE

def calculateJPLOEs(jploeDict, deltaCenturies):
    a = calculateJPLOE(jploeDict["a0"], jploeDict["aDot"], deltaCenturies, cF=False) # au, au/cty
    ecc = calculateJPLOE(jploeDict["e0"], jploeDict["eDot"], deltaCenturies, cF=False) # unitless
    I = calculateJPLOE(jploeDict["I0"], jploeDict["IDot"], deltaCenturies, cF=False)
    LOP = calculateJPLOE(jploeDict["LOP0"], jploeDict["LOPDot"], deltaCenturies, cF=False) # deg, deg/cty
    LAN = jploeDict["LAN0"] + deltaCenturies * jploeDict["LANDot"] # deg. deg/cty
    if "cF" in jploeDict:
        cF = jploeDict["cF"]
        L = calculateJPLOE(jploeDict["L0"], jploeDict["LDot"], deltaCenturies, cF)
    else:
        L = calculateJPLOE(jploeDict["L0"], jploeDict["LDot"], deltaCenturies) # deg, deg/cty

    I = np.deg2rad(I)
    LOP = np.deg2rad(LOP)
    LAN = np.deg2rad(LAN)
    return a, ecc, I, L, LOP, LAN

def jploeToPerifocal(a, ecc, I, L, LDot, LOP, LAN):
    eccSquare = ecc**2
    argOfPeri = LOP - LAN
    M = L - LOP
    M = mappingUtils.radModPlusMinusPi(M)
    E = solveKepler(M, ecc)
    EDot = LDot/(1.0 - ecc * np.cos(E))
    eccSqrtTerm = np.sqrt(1.0 - eccSquare)

    xPeri = a * (np.cos(E) - ecc) # Object's coordinate - x axis aligned from focus to periapsis
    vxPeri = -a * np.sin(E) * EDot
    yPeri = a * eccSqrtTerm * np.sin(E)
    vyPeri = a * np.cos(E) * EDot * eccSqrtTerm

    print("vxPeri:", vxPeri)
    print("vyPeri:", vyPeri)

    periPos = np.array([xPeri, yPeri])
    periVel = np.array([vxPeri, vyPeri])

    radius = np.sqrt(xPeri**2.0 + yPeri**2.0)
    print("R:", radius)
    print()
    return periPos, periVel

def jploeToCartesian2D(a, ecc, L, LDot, LOP, LAN):
    I = 0.0
    periPos, periVel = jploeToPerifocal(a, ecc, I, L, LDot, LOP, LAN)
    argOfPeri = LOP - LAN
    x = periPos[0] * np.cos(argOfPeri) - periPos[1] * np.sin(argOfPeri)
    vx = periVel[0] * np.cos(argOfPeri) - periVel[1] * np.sin(argOfPeri)
    y = periPos[0] * np.sin(argOfPeri) + periPos[1] * np.cos(argOfPeri)
    vy = periVel[0] * np.sin(argOfPeri) + periVel[1] * np.cos(argOfPeri)
    pos = np.array([x, y])
    vel = np.array([vx, vy])
    return pos, vel

def jploeToCartesian3D(a, ecc, EA, AOP, I, LAN, LDot): #inputs in /cty, rad and rad/cty, output in AU and AU per century
    ecc2 = ecc**2
    P = a * (np.cos(EA) - ecc)
    Q = a * np.sin(EA) * np.sqrt(1 - ecc2)

    vP = -a * np.sin(EA) * LDot / (1 - ecc * np.cos(EA))
    vQ = a * np.cos(EA) * np.sqrt(1 - ecc2) * LDot / (1 - ecc * np.cos(EA))

    x = np.cos(AOP) * P - np.sin(AOP) * Q
    vX = np.cos(AOP) * vP - np.sin(AOP) * vQ

    y = np.sin(AOP) * P + np.cos(AOP) * Q
    vY = np.sin(AOP) * vP + np.cos(AOP) * vQ

    #Assume Inclination = 0.0 - two dimensional
    z = np.sin(I) * x
    vZ = np.sin(I) * vX

    x = np.cos(I) * x
    vX = np.cos(I) * vX

    xT = x
    vXT = vX

    x = np.cos(LAN) * xT - np.sin(LAN) * y
    vX = np.cos(LAN) * vXT - np.sin(LAN) * vY

    y = np.sin(LAN) * xT + np.cos(LAN) * y
    vY = np.sin(LAN) * vXT + np.cos(LAN) * vY

    position = np.array([x, y, z])
    velocity = np.array([vX, vY, vZ])

    return position, velocity

def getSVFromJPLOE(jploeDict, initTime, parent):
    a, ecc, I, L, LOP, LAN = getAllCurrentOEs(jploeDict, initTime)
    AOP = LOP - LAN
    AOPRad = np.deg2rad(AOP)
    M = L - LOP
    M = degModuloHalfTopHalfBot(M)
    MRad = np.deg2rad(M)
    IRad = np.deg2rad(I)
    LANRad = np.deg2rad(LAN)
    LDot = jploeDict["LDot"]
    LDotRad = np.deg2rad(LDot)
    E = solveKepler(MRad, ecc)

    pos, vel = perifocalToCartesian(a, ecc, E, AOPRad, IRad, LANRad, LDotRad)
    pos = pos * au
    vel = vel * au / (100 * 365.25 * 24 * 60 * 60) #Convert from au / cty to m / s
    print("Parent Pos:", parent.pos)
    print("Parent Vel:", parent.vel)
    if not (parent == None):
        pos += parent.pos
        vel += parent.vel
    return pos, vel

def getOEFromSV(position, velocity, parentObj):
    eps = 1e-6
    if (parentObj == None):
        return False #if not orbiting anything, doesn't have parent object
    else:
        relPos = self.pos - parent.pos
        radius = np.linalg.norm(relPos)
        relVel = self.vel - parent.vel
        vSquared = relvel.dot(relVel)
        u = parent.gravParam
        specMechEng = (vSquared/2) - (uOverR)

        return calculateOE(relPos, relVel, u, specMechEng)

def calculateEnergies(radVec, velVec, gravParam):
    h = np.cross(radVec, velVec)
    vSquared = velVec.dot(velVec)
    r = np.linalg.norm(radVec)
    specMechEng = vSquared/2.0 - gravParam/r
    return h, specMechEng

def calculateOE2D(radVec, velVec, gravParam, specMechEng):
    radius = np.linalg.norm(radVec)
    muOverR = gravParam/radius
    rDotV = radVec.dot(velVec)
    vSquared = velVec.dot(velVec)
    nMag = 0.0 # n = k cross h = [0.0, 0.0, 1.0] cross [0.0, 0.0, hMag]

    eccVec = ((vSquared - muOverR)*radVec - rDotV * velVec)/gravParam
    e = np.linalg.norm(eccVec)

    eps = 1e-6
    if (np.abs(e - 1.0) <= eps): #Eccentricity close to or greater than 1 -  Parabolic
        #p = (h.dot(h))/gravParam
        a = np.inf
    else:
        a = -(gravParam / (2 * specMechEng))
        #p = a * (1 - eccMag**2)

    if (e == 0):
        argOfPeri = 0.0
        trueAnomaly = np.arctan2(radVec[1], radVec[0])
    else:
        argOfPeri = np.arctan2(eccVec[1], eccVec[0]) #
        trueAnomaly = np.arccos(eccVec.dot(radVec)/(e * radius))

    if(rDotV < 0.0):
        trueAnomaly = (2 * np.pi - trueAnomaly)
    return a, e, argOfPeri, trueAnomaly #I and LAN are both taken to be 0 degrees

def calculateOE3D(radVec, velVec, gravParam, specMechEng):
    h = np.cross(radVec, velVec)
    #print("h:", h)
    hMag = np.linalg.norm(h)
    zVec = np.array([0.0, 0.0, 1.0])
    nVec = np.cross(zVec, h)
    nMag = np.linalg.norm(nVec)
    rDotV = radVec.dot(velVec)
    vSquared = velVec.dot(velVec)
    radius = np.linalg.norm(radVec)
    muOverR = gravParam/radius
    eccVec = ((vSquared - muOverR)*radVec - (rDotV)*velVec)/gravParam
    eccMag = np.linalg.norm(eccVec)

    eps = 1e-6
    if (np.abs(eccMag - 1.0) <= eps): # Circular
        p = (h.dot(h))/gravParam
        a = np.inf
    else:
        a = -(gravParam / (2 * specMechEng))
        p = a * (1 - eccMag**2)

    hK = h[2]
    eK = eccVec[2]

    eps = 1e-6
    if (nMag < eps):
        Omega = 0.0
    else:
        nI = nVec[0]
        nJ = nVec[1]
        Omega = np.arccos(nI/nMag)
        if (nJ < 0.0):
            Omega = (2*np.pi - Omega)

    i = np.arccos(hK/hMag)

    #argOfPeri = np.arccos(nVec.dot(eccVec)/(nMag * eccMag))
    argOfPeri = 0.0
    trueAnomaly = np.arccos(eccVec.dot(radVec)/(eccMag * radius))

    if (eK < 0.0):
        argOfPeri = (2*np.pi - argOfPeri)

    if (rDotV < 0.0):
        trueAnomaly = (2*np.pi - trueAnomaly)

    return a, eccMag, i, Omega, argOfPeri, trueAnomaly

def getROI(smallerMass, largerMass, a):
    powerFactor = (2.0/5.0)
    roiRatio = (smallerMass / largerMass) ** powerFactor
    roi = roiRatio * a
    return roi
