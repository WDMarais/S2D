import numpy as np
import math
import json
import sys
import datetime

from scipy.constants import G # Newton's gravitational constant
from scipy.constants import au

import timeUtils
import orbitUtils
import fileUtils
import sscUtils
import mappingUtils

from constants import IDEAL_PSOLAR
from miscConstants import DIMENSION_ERROR
from miscConstants import LAZY_CODER_ERROR

genPath = "general.json"
with open(genPath) as g:
    generalSettings = json.load(g)
g.close()

scenePath = generalSettings["sceneFile"]
with open(scenePath) as s:
    scene = json.load(s)
s.close()

dimensions = scene["numDimensions"]

orbits = {
            "InitOrbit":0,
            "PlanetOrbit":1,
            "SunOrbit":2,
            "TransferOrbit":3,
            "TargetOrbit":4,
            "ExitPrep":5,
            "EntryPrep":6
        }

class spaceBody(object):
    def __init__(self, name, physicsDict, initTime = None, parent = None):
        self.name = name
        self.mass = physicsDict["mass"]
        self.gravParam = self.mass * G
        self.orbitParent = parent #Assumption from settings
        if (initTime == None):
            initTime = timeUtils.returnJ2000()
        initDict = physicsDict["initialization"]
        ctg = initDict["category"]
        if (ctg == "CARTESIAN"):
            self.setStateFromSV(initDict["pos"], initDict["vel"], parent)
        elif (ctg == "JPLOE"):
            if(dimensions == 2):
                self.setStateFromJPLOE2D(initTime, initDict, parent)
            elif(dimensions == 3):
                print("init SpaceBody - Three dimensions")
                sys.exit(LAZY_CODER_ERROR)
        elif (ctg == "KOE2D"):
            if(dimensions == 2):
                self.setStateFromKOE2D(initDict)
            elif(dimensions == 3):
                print("init SpaceBody - Three dimensions")
                sys.exit(LAZY_CODER_ERROR)

        self.forces = np.zeros(dimensions)
        self.rot = 0

    def setStateFromSV(self, pos, vel, parent = None):
        posDim = len(pos)
        velDim = len(vel)
        if not(posDim == dimensions):
            print("Position Dimensioning Error:")
            print("Expected", dimensions)
            print("Received", posDim)
            sys.exit(DIMENSION_ERROR)
        elif not(velDim == dimensions):
            print("Velocity Dimensioning Error:")
            print("Expected", dimensions)
            print("Received", velDim)
            sys.exit(DIMENSION_ERROR)
        else:
            npPos = np.array(pos)
            npVel = np.array(vel)
            if not (parent == None):
                npPos = npPos + np.array(parent.pos)
                npVel = npVel + np.array(parent.vel)
            self.pos = npPos
            self.vel = npVel

    def setStateFromJPLOE2D(self, initTime, jploeDict, parent = None):
        deltaCenturies = timeUtils.centuriesSinceJ2000(initTime)
        a, ecc, I, L, LOP, LAN, = orbitUtils.calculateJPLOEs(jploeDict, deltaCenturies)
        LDot = np.deg2rad(jploeDict["LDot"])
        print("spaceBodies setStateFromJPLOE2D", self.name)
        pos, vel = orbitUtils.jploeToCartesian2D(a, ecc, L, LDot, LOP, LAN)

        distScale = au #Given in au, multiply by au unit - m/au
        timeScale = timeUtils.secondsToCenturies(1.0)
        pos *= distScale
        vel *= (distScale*timeScale)

        self.pos = pos
        self.vel = vel

        if not (parent == None):
            self.pos += parent.pos
            self.vel += parent.vel

    def updateStateFromKOE2D(self, rotDir):
        p = self.a * (1 - (self.e)**2.0)
        nu = self.trueAnomaly
        gravParam = self.gravParam + self.orbitParent.gravParam

        posBase = p / (1 + self.e * np.cos(nu))
        posFracts = np.array([np.cos(nu), np.sin(nu)])
        posPQ = posBase * posFracts

        velBase = np.sqrt(gravParam/p)
        velFracts = np.array([-np.sin(nu), (self.e + np.cos(nu))])
        velPQ = velBase * velFracts

        w = self.argOfPeri
        pos = mappingUtils.rotate2D(posPQ, w)
        vel = mappingUtils.rotate2D(velPQ, w)

        if (rotDir == "cw"):
            vel *= -1.0

        if (self.orbitParent != None):
            pos += self.orbitParent.pos
            vel += self.orbitParent.vel

        self.pos = pos
        self.vel = vel

    def setStateFromKOE2D(self, oeDict):
        self.a = oeDict["a"] * au
        self.e = oeDict["e"]
        self.argOfPeri = oeDict["argOfPeri"] * (np.pi/180)
        self.trueAnomaly = oeDict["trueAnomaly"] * (np.pi/180)
        rotDir = oeDict["rotDir"]
        self.updateStateFromKOE2D(rotDir)

    def setTarget(self, target):
        self.orbitTarget = target
        self.orbitType = orbits["InitOrbit"] #Assumption, updated in next steps

    def calculateFg(self, other):
        rVec = other.pos - self.pos
        rSquared = np.dot(rVec, rVec)
        rMag = np.sqrt(rSquared)
        fDir = rVec / rMag
        f = self.gravParam * (other.mass / rSquared) * fDir
        return f

    def resetForces(self):
        self.forces = np.zeros(dimensions)

    def setOrbitParent(self, newOrbitParent):
        self.orbitParent = newOrbitParent
        if (self.orbitParent == None):
            print(self.name, "has no parent")
        else:
            print(self.name, "new parent:", self.orbitParent.name)

    def computeGravity(self, others):   # Computes and assigns net gravitational
        fMax = 0.0                      # forces on each body in body list.
        for other in others:
            f = self.calculateFg(other)
            self.forces += f

    def updateSV(self, dT):
        self.acc = self.forces / self.mass
        self.vel += dT * self.acc
        self.pos += dT * self.vel

    def updateOEs2D(self):
        op = self.orbitParent
        if (op == None):
            self.a = 0.0
            self.e = 0.0
            self.argOfPeri = 0.0
            self.trueAnomaly = 0.0
            self.h = 0.0
        else:
            r = self.pos - op.pos
            v = self.vel - op.vel
            gravParam = self.gravParam + op.gravParam
            self.h, self.specMechEng = orbitUtils.calculateEnergies(r, v, gravParam)
            a, e, argOfPeri, trueAnomaly = orbitUtils.calculateOE2D(r, v, gravParam, self.specMechEng)
            self.a = a
            self.e = e
            self.argOfPeri = argOfPeri
            self.trueAnomaly = trueAnomaly

    def updateAuxVars(self):
        self.p = self.a * (1 - self.e**2.0)
        if (self.orbitParent == None):
            self.meanMotion = 0.0
        else:
            self.meanMotion

    def getOEs2D(self):
        return self.a, self.e, self.argOfPeri, self.trueAnomaly

    def updateROI(self):
        if (self.orbitParent == None):
            self.roi = np.inf
        else:
            self.roi = self.a * ((self.mass / self.orbitParent.mass)**(2.0/5.0))

    def checkPossibleParent(self, other):
        if (self.mass > other.mass):
            isPossible = False
        else:
            rVec = self.pos - other.pos
            rMag = np.linalg.norm(rVec)
            bodyOutsideRoi = (rMag > other.roi)
            if bodyOutsideRoi:
                isPossible = False
            else:
                isPossible = True
        return isPossible

    def updateParent(self, sun, others):
        if (self != sun):
            parent = self.orbitParent
            parentRoiPlus = (1.05 * parent.roi) #Extend sphere slightly for "hysteresis curve" behaviour
            rVec = self.pos - parent.pos
            rMag = np.linalg.norm(rVec)
            if (rMag > parentRoiPlus): #Have moved out of parent SOI
                for other in others:
                    if self.checkPossibleParent(other):
                        self.setOrbitParent(other)
            else:#May have entered a new SOI smaller than parent SOI
                for other in others:
                    if (other.roi < parent.roi):
                        if self.checkPossibleParent(other):
                            self.setOrbitParent(other)

    def updateTrackingVariables(self, sun, others):
        self.updateOEs2D()
        self.updateROI()
        self.updateParent(sun, others)
        #print("updateTrackingVariables")

    def initArrays(self, numFrames):
        self.posArray = np.zeros((numFrames, dimensions))
        self.velArray = np.zeros((numFrames, dimensions))
        self.aArray = np.zeros((numFrames))
        self.eArray = np.zeros((numFrames))
        self.argOfPeriArray = np.zeros((numFrames))
        self.trueAnomalyArray = np.zeros((numFrames))
        self.hArray = np.zeros((numFrames))
        if (dimensions == 2):
            self.rotArray = np.zeros((numFrames))
        elif (dimensions == 3):
            self.rotArray = np.zeros((numFrames, numDimensions))
            self.iArray = np.zeros((numFrames))
            self.LANArray = np.zeros((numFrames))
        else:
            print("Wrong Dimensions in spaceBodies initArrays")
            sys.exit(INVALID_PARAMETER_ERROR)

    def storeStates(self, timeIndex):
        self.posArray[timeIndex] = self.pos
        self.velArray[timeIndex] = self.vel
        self.rotArray[timeIndex] = self.rot
        self.hArray[timeIndex] = self.h
        self.aArray[timeIndex] = self.a
        self.eArray[timeIndex] = self.e
        self.argOfPeriArray[timeIndex] = self.argOfPeri
        self.trueAnomalyArray[timeIndex] = self.trueAnomaly
        if (dimensions == 3):
            self.iArray[timeIndex] = self.i
            self.LANArray[timeIndex] = self.LAN

    def storeToFile(self, fileDir, timeArray):
        fullDir = fileDir + '/' + self.name + '/'
        fileUtils.storeJoinedArrays(fullDir, "Pos.dat", timeArray, self.posArray)
        fileUtils.storeJoinedArrays(fullDir, "Vel.dat", timeArray, self.velArray)
        fileUtils.storeJoinedArrays(fullDir, "Rot.dat", timeArray, self.rotArray)
        fileUtils.storeJoinedArrays(fullDir, "h.dat", timeArray, self.hArray)
        oeArray = np.column_stack((self.aArray, self.eArray, self.argOfPeriArray, self.trueAnomalyArray))
        if (dimensions == 3):
            oeArray = np.column_stack((oeArray, self.iArray, self.LANArray))
        fileUtils.storeJoinedArrays(fullDir, "OE.dat", timeArray, oeArray)
        print(self.name, " data saved.")

    def alignToPos2D(self, location):
        relPos = location - self.pos
        self.sailDir = relPos / np.linalg.norm(relPos)
        self.n = -self.sailDir
        self.rot = mappingUtils.dirToRadians2D(self.sailDir)
        # transverse direction points to sail "right" (not "left"),
        # relative to sail norm and orbit plane
        tRadians = self.rot + (np.pi/2.0)
        self.t = mappingUtils.radiansToDir2D(tRadians)

    def alignToRad2D(self, radians):
        self.rot = radians
        self.sailDir = mappingUtils.radiansToDir2D(radians)
        self.n = -self.sailDir
        # transverse direction points to sail "right" (not "left"),
        # relative to sail norm and orbit plane
        tRadians = self.rot + (np.pi/2.0)
        self.t = mappingUtils.radiansToDir2D(tRadians)

    def alignToDir2D(self, newDir):
        self.sailDir = newDir
        self.n = -newDir
        self.rot = mappingUtils.dirToRadians2D(newDir)
        # transverse direction points to sail "right" (not "left"),
        # relative to sail norm and orbit plane
        tRadians = self.rot + (np.pi/2.0)
        self.t = mappingUtils.radiansToDir2D(tRadians)

class forceSatellite(spaceBody):
    # Gets to utilise arbitrary force to get where it needs to go
    # 2 approaches to changing x:
    #   Set direction such that x change is maximised
    #   Set direction such that changes other than x are minimised
    def setMaxForce(force):
        self.maxForce = force

    def setDesiredOrbit2D(a, ecc, argOfPeri, orbitParent):
        self.aTarget = a
        self.eTarget = ecc
        self.argOfPeriTarget = argOfPeri
        self.orbitParentTarget = orbitParent

class solarSailSatellite(spaceBody):
    def __init__(self, name, physicsDict, initTime = None, parent = None):
        super().__init__(name, physicsDict, initTime, parent)
        params= physicsDict["sailParams"]
        self.reflectivity = params["reflectivity"]
        self.specular = params["specular"]
        self.epsF = params["epsF"]
        self.epsB = params["epsB"]
        self.Bf = params["Bf"]
        self.Bb = params["Bb"]
        self.a = 1 - self.reflectivity

    def setSailConditions(self, area, P_solar = IDEAL_PSOLAR, rho_s = 1.0, sailDir = None):
        self.sailArea = area
        self.reflectivity = rho_s
        self.P_solar = P_solar
        self.PA = self.sailArea * self.P_solar
        if (sailDir == None):
            self.sailDir = np.array([1.0, 0.0]) # pointing right
        else:
            sailDir = sailDir / np.linalg.norm(sailDir)
            self.sailDir = np.array(sailDir)
        self.n = -self.sailDir

    def computeIdealFs2D(self, sun):
        rayVec = self.pos - sun.pos
        rayDir = rayVec / np.linalg.norm(rayVec)
        rayDotNormal = rayDir.dot(self.n)
        if((rayDotNormal - 1.0) < 1e-6):
            self.pitchAngle = 0.0 #rounding errors, r dot n approx 1, arccos approx 0
        elif((rayDotNormal + 1.0) > 1e-6):
            self.pitchAngle = np.pi #rounding errors, r dot n approx -1, arccos approx PI rad
        else:
            self.pitchAngle = np.arccos(rayDotNormal)

        if (self.pitchAngle >= np.pi/2):
            self.solarForce = 0.0
        else:
            exposureScale = ((rayDir.dot(self.n))**2.0) * self.n
            self.solarForce = 2 * self.PA * exposureScale

        self.forces += self.solarForce

    def computeOpticFs2D(self, sun):
        rayVec = self.pos - sun.pos
        rayDir = rayVec / np.linalg.norm(rayVec)
        rayDotNormal = rayDir.dot(self.n)
        self.pitchAngle = np.arccos(rayDotNormal)
        print("spaceBodies computeOpticFs2D NOT YET COMPLETED")
        system.exit(LAZY_CODER_ERROR)
        if (self.pitchAngle >= np.pi/2):
            self.solarForce = 0.0
        else:
            pitch = self.pitchAngle
            cosOfPitch = np.cos(pA)
            sinOfPitch = np.sin(pA)
            u = (sinOfPitch * self.t) + (cosOfPitch * self.n)
            s = (sinOfPitch * self.t) - (cosOfPitch * self.n)

            fan = self.PA * (cosOfPitch ** 2.0) * self.n
            fat = self.PA * (cosOfPitch * sinOfPitch) * self.t
            fa = fan + fat

            frs = -1.0 # Complete
            fru = self.Bf * self.reflectivity * ()

    '''
    def widenOrbitAround(self, orbitBody):
        r = orbitBody.pos - self.pos
        rU = np.linalg.norm(r)
        n = self.sailDir
        p = rU.cross(n)
        q = np.linalg.norm(self.vel)

    def narrowOrbitAround(self, orbitBody):
        r = orbitBody.pos - self.pos
        rU = np.linalg.norm(r)
        n = self.sailDir
        p = rU.cross(n)
        q = np.linalg.norm(self.vel)
    '''

    def fDirLargestEccChange2D(self, sun):
        lineFromSun = self.pos - sun.pos

    def turnFromSun(self, sun, radAway):
        lineToSun = sun.pos - self.pos
        lsDir = lineToSun / np.linalg.norm(lineToSun)
        #self.alignToPos2D(lsDir)
        lsRad = mappingUtils.dirToRadians2D(lsDir)
        lsRad += (np.pi/6)
        lsRad = mappingUtils.radModPlusMinusPi(lsRad)
        self.alignToRad2D(lsRad)

    def raiseOrbitEnergy2D(self, sun):
        if (self.orbitParent == None):
            pass
        else:
            relVel = self.vel - self.orbitParent.vel
            rayLine = self.pos - sun.pos
            rayDir = rayLine / np.linalg.norm(rayLine)
            fDir = relVel / np.linalg.norm(relVel) #Max force added to v direction
            coneRad = sscUtils.coneRadForMaxF2D(rayDir, fDir)
            normDir = np.cross(rayDir, fDir)
            if(normDir < 0.0): # Pointing "down"
                coneRad = -coneRad #opposite direction
            fDirNew = sscUtils.coneRadToDir2D(rayDir, coneRad)
            sailDirNew = -fDirNew
            self.alignToDir2D(sailDirNew)
            '''
            self.alignToDir2D(fDir)
            '''
    def dropOrbitEnergy2D(self, sun):
        if (self.orbitParent == None):
            pass
        else:
            relVel = self.vel - sun.vel
            rayLine = self.pos - sun.pos
            rayDir = rayLine / np.linalg.norm(rayLine)
            fDir = -relVel / np.linalg.norm(relVel) #Max force added to v direction
            coneRad = sscUtils.coneRadForMaxF2D(rayDir, fDir)
            normDir = np.cross(rayDir, fDir)
            if(normDir < 0.0): # Pointing "down"
                coneRad = -coneRad #opposite direction
            fDirNew = sscUtils.coneRadToDir2D(rayDir, coneRad)
            sailDirNew = -fDirNew
            self.alignToDir2D(sailDirNew)
            '''
            self.alignToDir2D(fDir)
            '''

    def controllerStep(self, t, dT, sun):
        #self.turnFromSun(sun, (np.pi/6))
        self.raiseOrbitEnergy2D(sun)
        #self.dropOrbitEnergy2D(sun)
        '''
        if (t % 100 == 0):
            relV = self.vel - self.orbitParent.vel
            relVDir = relV / np.linalg.norm(relV)
            #print(mappingUtils.dirToRadians2D(relVDir)*180/np.pi)
        '''
        self.computeIdealFs2D(sun)

class thrustSatellite(spaceBody):
    def setThrustVector(self, thrustVector=np.zeros(dimensions)):
        self.thrustVector = thrustVector

    def controllerStep(self, t, dT):
        pass
