from spaceBodies import spaceBody
from spaceBodies import thrustSatellite
from spaceBodies import solarSailSatellite
import numpy as np
from math import radians

from miscConstants import INVALID_PARAMETER_ERROR

def initBody(bodyDict, initTime, parent = None):
    name = bodyDict["name"]
    p = bodyDict["physics"]
    c = bodyDict["category"]

    if (c =="GENERAL"):
        body = spaceBody(name, p, initTime, parent)
    elif (c =="THRUSTSAT"):
        body = thrustSatellite(name, p, initTime, parent)
    elif (c =="SOLSAT"):
        body = solarSailSatellite(name, p, initTime, parent)
        area = p["sailArea"]
        body.setSailConditions(area) #Other parameters kept at defaults
    else:
        print("simUtils initBody Error - Category", c, "Invalid")
        sys.exit(1)
    return body

def initBodiesRecursive(bodiesDict, initTime, parent = None):
    bodies = []
    nameMap = {}
    for b in bodiesDict:
        body = initBody(b, initTime, parent)
        #body.logState()
        bodies = bodies + [body]
        nameMap.update({body.name:body})
        if "children" in b:
            children, cNameMap = initBodiesRecursive(b["children"], initTime, body)
            bodies = bodies + children
            nameMap.update(cNameMap)
    return bodies, nameMap

def flattenDict(bodiesDict):
    flattenedDict = []
    for b in bodiesDict:
        if "children" in b:
            children = b["children"]
            del b["children"]
            flattenedDict = flattenedDict + [b]
            childrenDict = flattenDict(children)
            flattenedDict = flattenedDict + childrenDict
        else:
            flattenedDict = flattenedDict + [b]
    return flattenedDict

def getHeaviestBody(bodies):
    maxMass = -np.inf
    result = None
    for body in bodies:
        if (body.mass > maxMass):
            maxMass = body.mass
            result = body

    if (result == None):
        print("simUtils getHeaviestBody Error: no bodies with appropriate mass")
        sys.exit(INVALID_PARAMETER_ERROR)

    return result
