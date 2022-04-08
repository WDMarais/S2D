import bpy, bmesh
import json
from scipy.constants import au
import numpy as np
import mathutils
from math import radians

#####################
# Blender Utilities #
#####################
O = bpy.ops
C = bpy.context
D = bpy.data

genPath = "general.json"
colorsPath = "colors.json"

with open(colorsPath) as colors:
    c = json.load(colors)
colors.close()

with open(genPath) as g:
    generalSettings = json.load(g)
g.close()

scenePath = generalSettings["sceneFile"]

with open(scenePath) as s:
    scene = json.load(s)
s.close()

sceneName = scene["name"]
bodiesDict = scene["bodies"]
numSteps = scene["numSteps"]

trailDuration = scene["trailDuration"]
numParticles = scene["numParticles"]

def scaleDistFromSI(newUnit):
    if (newUnit == 'm'):
        scaleFactor = 1.0
    elif (newUnit == 'au'):
        scaleFactor = (1.0 / au)
    elif (newUnit == 'km'):
        scaleFactor = (1e-3)
    else:
        print("render.py scaleDistFromSI error - invalid parameter")
        print(newUnit)
        sys.exit(INVALID_PARAMETER_ERROR)
    return scaleFactor

def flattenRecursiveBodiesDict(bodiesDict):
    print("Loading body JSONs")
    flattenedDict = []
    for b in bodiesDict:
        if "children" in b:
            children = b["children"]
            del b["children"]
            flattenedDict = flattenedDict + [b]
            childrenDict = flattenRecursiveBodiesDict(children)
            flattenedDict = flattenedDict + childrenDict
        else:
            flattenedDict = flattenedDict + [b]
    return flattenedDict

def setMaterial(obj, mat):
    obj.active_material = mat

def createMaterial(name, color, whichType):
    mat = D.materials.get(name)
    if mat is None:
        mat = bpy.data.materials.new(name)
    if len(color) == 3:
        cList = list(color)
        cList.append(1.0) ## add alpha channel - no transparency, fully opaque
        color = tuple(cList)
    mat.diffuse_color = color
    mat.specular_intensity = 1.0
    mat.type = whichType
    mat.emit = 0.5
    return mat

def clearAllObjects():
    for o in bpy.data.objects:
        o.select_set(True)
    bpy.ops.object.delete(use_global=False)

    for m in bpy.data.materials:
        bpy.data.materials.remove(m)

def addVertexObj(name, location):
    O.mesh.primitive_plane_add(location = location)
    opObj = O.object
    opObj.mode_set(mode="EDIT")
    O.mesh.select_all(action="DESELECT")
    opObj.mode_set(mode="OBJECT")

    vObj = C.object
    vObj.name = name

    verts = vObj.data.vertices
    for i in range(1, 4):
        verts[i].select_set(True)
    opObj.mode_set(mode="EDIT")
    O.mesh.delete(type="VERT")
    opObj.mode_set(mode="OBJECT")
    opObj.origin_set(type="GEOMETRY_ORIGIN")

    return vObj

def addParticles(obj, start, end, lifetime):
    particleSysName = obj.name + "Particles"
    C.scene.objects.active = obj
    if (particleSysName not in obj):
        O.object.particle_system_add()
        C.object.particle_systems[0].name = particleSysName

    particleSys = obj.particle_systems[particleSysName]
    settings = particleSys.settings
    settings.emit_from = "VERT"
    settings.effector_weights.gravity = 0
    settings.frame_start = start
    settings.frame_end = end
    settings.lifetime = lifetime
    settings.count = numParticles

def addTrail(parentObj, trailWidth):
    start = 1
    numSteps = scene["numSteps"]
    location = (0.0, 0.0, 0.0)
    trailName = parentObj.name + "Trail"

    vert = addVertexObj(trailName, location)
    vert.parent = parentObj
    parentColor = parentObj.active_material.diffuse_color
    vertMat = createMaterial(trailName + "Material", parentColor, "HALO")
    vertMat.halo.size = trailWidth
    setMaterial(vert, vertMat)
    addParticles(vert, start, numSteps, trailDuration) #Rendering all frames, so the duration and the end are the same

def addSolarSail(name, size, location = (5, 5, 5)):
    sailDisplacement = 1.5 * size * np.array([1.0, 0.0, 0.0])
    #sailThickness = tuple(0.2 * size * np.array([0.0, 0.0, 1]))
    O.mesh.primitive_cylinder_add(radius=2*size, depth=size)
    A = C.active_object
    A.scale = [1.0, 1.0, 0.2]
    O.transform.rotate(value=(-np.pi/2), axis=(0.0, 1.0, 0.0), constraint_axis=(False, True, False), constraint_orientation='GLOBAL', mirror=False, proportional='DISABLED', proportional_edit_falloff='SMOOTH', proportional_size=1)
    A.location = location + sailDisplacement
    O.mesh.primitive_uv_sphere_add(scale=(size, size, size), location=location)
    B = C.active_object
    O.object.shade_smooth()
    B.name = name

    A.select_set(True)
    B.select_set(True)
    O.object.join()

def addActiveBody(name, location, graphics, blenderUnits):
    r = graphics["representation"]
    s = graphics["size"] * blenderUnits
    c = graphics["color"]

    with open("colors.json") as cJSON:
        cList = json.load(cJSON)
    cJSON.close()

    color = cList[c]

    if (r == "SPHERE"):
        O.mesh.primitive_uv_sphere_add(scale=(s, s, s), location=location)
        O.object.shade_smooth()
    elif (r == "SOLSAT"):
        addSolarSail(name, s, location)

    body = bpy.context.active_object
    body.name = name
    matName = name + "Material"
    mat = createMaterial(matName, color, "SURFACE")
    setMaterial(body, mat)

    leaveTrail = (graphics["leaveTrail"] == "True")
    print(name, leaveTrail)
    if leaveTrail:
        trailWidth = graphics["trailWidth"]
        addTrail(body, trailWidth)

def addCamera(name, loc, rot, parent=None):
    bpy.ops.object.camera_add(location = loc, rotation = rot)
    body = bpy.context.active_object
    body.name = name
    body.location = loc
    body.rotation_euler = rot
    if not(parent == None):
        body.parent = parent

def initializeScene(firstFrame, lastFrame):
    clearAllObjects()
    bpy.context.scene.frame_start = firstFrame
    bpy.context.scene.frame_end = lastFrame

def addKeyFrame(objectName, objectLocation):
    o = bpy.data.objects[objectName]
    o.location = objectLocation
    o.keyframe_insert(data_path="location", index=-1)

def setRenderSettings(settingDict):
    bpy.context.scene.render.resolution_x = settingDict["resolution_x"]
    bpy.context.scene.render.resolution_y = settingDict["resolution_y"]
    bpy.context.scene.render.resolution_percentage = settingDict["resolution_percentage"]
    bpy.context.scene.render.ffmpeg.codec = settingDict["file_codec"]
    extension = ".mp4"
    bpy.context.scene.render.filepath = "//movies/" + sceneName + extension
    bpy.context.scene.render.ffmpeg.format = settingDict["file_format"]
    bpy.context.scene.render.fps = settingDict["frames_per_second"]
    print("Rendering settings setup")

#EXECUTED CODE
setRenderSettings(generalSettings)
dimensions = scene["numDimensions"]
flatBD = flattenRecursiveBodiesDict(bodiesDict)
numBodies = len(flatBD)
numPosCoords = dimensions+1
if (dimensions == 2):
    numRotCoords = 2
    raShape = (numBodies, numSteps, numRotCoords)
elif (dimensions == 3):
    numRotCoords = 4
    raShape = (numBodies, numSteps, numRotCoords)

paShape = (numBodies, numSteps, numPosCoords)
posArray = np.zeros(paShape)
rotArray = np.zeros(raShape)

sceneDataPath = generalSettings["storeParentDir"] + "/" + sceneName
for b, body in enumerate(flatBD):
    name = body["name"]
    posFilePath = f"{sceneDataPath}/{name}/Pos.dat"
    rotFilePath = f"{sceneDataPath}/{name}/Rot.dat"
    posArray[b] = np.genfromtxt(posFilePath, delimiter=',')
    rotArray[b] = np.genfromtxt(rotFilePath, delimiter=',')
    print("Loaded", name, "data")

viewCenter = scene["viewCenter"]
if not(viewCenter == 'Origin'):
    centerPath = f"{sceneDataPath}/{viewCenter}/Pos.dat"
    centerPosArray = np.genfromtxt(centerPath, delimiter=',')
    for b, body in enumerate(flatBD):
        posArray[b] -= centerPosArray

distStoreUnit = scene["distStoreUnit"]
distRenderUnit = scene["distRenderUnit"]
#baseUnitScale = scene["baseUnitScale"] #What is to be used as the distance for a single unit
blenderScale = scene["blenderScale"] #How many blender units is to correspond to a base unit

if (distStoreUnit != distRenderUnit):
    storeScale = scaleDistFromSI(distStoreUnit)
    renderScale = scaleDistFromSI(distRenderUnit)
    scaleChange = renderScale / storeScale
    posArray *= scaleChange

posArray *= blenderScale
#addBodies(flatBD, leaveTrails, posArray)

initializeScene(1, numSteps)
for b, body in enumerate(flatBD):
    name = body["name"]
    position = posArray[b][0]
    if (len(position) == 2):
        position = np.append(position, [0.0])
    g = body["graphics"]
    addActiveBody(name, position, g, blenderScale)
    print("Added", name, "render object")

bpy.context.scene.frame_current = 1
bpy.ops.object.select_all(action='TOGGLE')

def locKeyFrame(obj, objectLocation):
    x = objectLocation[1]
    y = objectLocation[2]
    z = 0.0
    obj.location = (x,y,z)
    obj.keyframe_insert(data_path="location")

def rotKeyFrame(obj, objectRotation):
    obj.rotation_euler[2] = objectRotation[1]
    obj.keyframe_insert(data_path="rotation_euler")

def noteLargeRadChange(oldRad, newRad):
    radChange = np.abs(oldRad - newRad)
    doNote = (radChange > (0.95 * np.pi))
    return doNote

def keyFrameBodies(frame):
    C.scene.frame_current = frame
    for b, body in enumerate(flatBD):
        name = body["name"]
        o = bpy.data.objects[name]
        pos = posArray[b][frame]
        rot = rotArray[b][frame]
        locKeyFrame(o, pos)
        rotKeyFrame(o, rot)
        if (frame < lastKeyFrame):
            rotNext = rotArray[b][frame+stepsPerKeyFrame]
            preventGimbal = noteLargeRadChange(rot[1], rotNext[1])
            if preventGimbal:
                for i in range(1, stepsPerKeyFrame):
                    firstFrame = frame + i
                    nextFrame = firstFrame + 1
                    firstRot = rotArray[b][firstFrame]
                    nextRot = rotArray[b][nextFrame]
                    gimbalJump = noteLargeRadChange(firstRot[1], nextRot[1])
                    if gimbalJump:
                        print("Patching")
                        C.scene.frame_current = firstFrame
                        rotKeyFrame(o, firstRot)
                        C.scene.frame_current = nextFrame
                        rotKeyFrame(o, nextRot)

def keyFrameBodies2(frame):
    frameIndex = frame // stepsPerKeyFrame
    C.scene.frame_current = frameIndex
    for b, body in enumerate(flatBD):
        name = body["name"]
        o = bpy.data.objects[name]
        pos = posArray[b][frame]
        rot = rotArray[b][frame]
        locKeyFrame(o, pos)
        rotKeyFrame(o, rot)

stepsPerKeyFrame = scene["stepsPerKeyFrame"]
lastKeyFrame = numSteps - stepsPerKeyFrame - 1
consoleUpdateFreq = scene["logFrequency"]

numFrames = (numSteps // stepsPerKeyFrame)
for n in range(numFrames):
    keyFrameBodies(n*stepsPerKeyFrame)

keyFrameBodies(numSteps - 1)
bpy.context.scene.frame_current = 1

name = "Tracker"
loc = np.array([0.0, 0.0, 0.05]) * blenderScale
rot = [0.0, 0.0, 0.0]
parent = D.objects["Earth"]
addCamera(name, loc, rot, parent)
Cam = D.objects[name]
Cam.data.clip_end = 20.0 * blenderScale

bpy.ops.object.lamp_add(type='SUN', view_align=False, location=(0, 0, 0))
