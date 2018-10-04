import simUtils
import timeUtils
import fileUtils
import mappingUtils
from spaceBodies import spaceBody
from spaceBodies import thrustSatellite
from spaceBodies import solarSailSatellite
import numpy as np
import json
from scipy.constants import au
import mathutils
from math import radians

from miscConstants import MAIN_JSON

def loadScene(genPath):
    sceneDict = {}
    bodyDict = {}
    genPath = "general.json"
    with open(genPath) as g:
        generalSettings = json.load(g)
    g.close()

    scenePath = generalSettings["sceneFile"]
    with open(scenePath) as s:
        scene = json.load(s)
    s.close()

    sceneDict.update({"storeParentDir":generalSettings["storeParentDir"]})
    sceneDict.update({"name":scene["name"]})
    sceneDict.update({"timeChange":scene["timeChange"]})
    sceneDict.update({"numSteps":scene["numSteps"]})
    sceneDict.update({"numDimensions":scene["numDimensions"]})
    sceneDict.update({"logFrequency":scene["logFrequency"]})
    sceneDict.update({"distStoreUnit":scene["distStoreUnit"]})
    sceneDict.update({"timeStoreUnit":scene["timeStoreUnit"]})
    sceneDict.update({"velStoreUnit":scene["velStoreUnit"]})
    sceneDict.update({"angleStoreUnit":scene["angleStoreUnit"]})
    initTime = timeUtils.dictToDateTime(scene["initTime"])
    bodies, bodyMap = simUtils.initBodiesRecursive(scene["bodies"], initTime)
    print("Loading body JSONs")
    flatBodyDict = simUtils.flattenDict(scene["bodies"])
    for f in flatBodyDict:
        p = f["physics"]
        if("orbitTarget" in p):
            body = bodyMap[f["name"]]
            oT = bodyMap[p["orbitTarget"]]
            body.setTarget(oT)

    for body in bodies:
        body.initArrays(scene["numSteps"])
    return sceneDict, bodies, bodyMap, initTime

def initBodiesStep(bodies, sun):
    for body in bodies:
        if (body != sun):
            body.orbitParent = sun
            body.updateOEs2D()
            body.updateROI()

    for b, body in enumerate(bodies):
        others = np.delete(bodies, b)
        body.updateTrackingVariables(sun, others)

def computationStep(bodies, t, dT, sun):
    for b, body in enumerate(bodies):
        others = np.delete(bodies, b)
        body.updateTrackingVariables(sun, others) # Orbit Elements, ROI, reference signals etc
        body.resetForces()
        body.computeGravity(others)
        bT = type(body)
        if (bT is thrustSatellite):
            body.controllerStep(t, dT)
        elif (bT is solarSailSatellite):
            body.controllerStep(t, dT, sun)

def updateStep(bodies, t, dT):
    for body in bodies:
        body.updateSV(dT)
        body.storeStates(t)

print("Startup")
sceneDict, bodies, bodyMap, initTime = loadScene(MAIN_JSON)
sun = simUtils.getHeaviestBody(bodies) #Assume sun is heaviest body

dT = sceneDict["timeChange"]
numSteps = sceneDict["numSteps"]
logFreq = sceneDict["logFrequency"]
tCount = logFreq
initBodiesStep(bodies, sun)
for t in range(numSteps):
    computationStep(bodies, t, dT, sun) # Calculate gravity, necessary feedback controllers etc
    updateStep(bodies, t, dT) # Based on previous step, update SVs, OEs and others
    if (tCount >= logFreq):
        print("Num Frames:", t)
        tCount = 1
    else:
        tCount += 1

unitScaleDict = {}

storeDistScale = mappingUtils.scaleDistFromSI(sceneDict["distStoreUnit"])
storeVelScale = mappingUtils.scaleVelFromSI(sceneDict["velStoreUnit"])
storeTimeScale = mappingUtils.scaleTimeFromSI(sceneDict["timeStoreUnit"])

timeArray = dT * np.linspace(0, (numSteps - 1), numSteps)
timeArray *= storeTimeScale

fileDir = sceneDict["storeParentDir"] + "/" + sceneDict["name"]
fileUtils.createDirs(bodies, fileDir)
for b in bodies:
    b.posArray *= storeDistScale
    b.aArray *= storeDistScale
    b.velArray *= storeVelScale
    b.storeToFile(fileDir, timeArray)
