from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
import random
import math

from ArcBall import *

name = "Dumb"

# by 400x400 square here, and frustrum later
# "400" pixels is one unit of virtual space (a 'virt').
# on this computer, I sit about 34 inches from the screen
# and 400 pixels is about 4.5 inches; so the frustum viewport
# should be 7.55 units away from the origin (my eye).
width = 400
height = 400
prec = 1
precinv = 1.0/prec

# radius in virtual units, for a sphere the size of the earth
# radius of earth is about 3,959 mi * 5280 ft/mi * 12 inch/ft *400pix/4.5in * 1 virt/400 pixels
#radius = 55742720 # virts
radius = 2000
virtsWidth = 1.0
virtsHeight = 1.0
virtsToEye = 7.55

displayListHandle = None

arcBall = None

# we want the initial distance to be such that, viewed as a circle on the screen initially
# the circle is centered on the window with radius z, and the window is 4z X 4z.
# Achieved via similar triangles: .5virtsWidth/virtsToEye == 2radius/initialDistance
initialDistance =  2 * radius*virtsToEye / virtsWidth

# distance from surface of sphere, can never get past surface of sphere
logDistanceE2 = int(math.log(initialDistance-radius)*100)

def init2():
    global arcBall
    ### generic ###
    glViewport(0,0,400,400)

    arcBall = ArcBall(400, 400)
    glutMouseFunc(lambda btn, state, x, y: arcBall.onClick(btn, state, x, y))
    glutMotionFunc(lambda x, y: arcBall.onMouseMove(x, y))

    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()
    #glOrtho(-50.0, 50.0, -50.0, 50.0, -50, 50.0)
    glFrustum(-virtsWidth*.5, virtsWidth*.5,
              virtsHeight*-.5, virtsHeight*.5,
              virtsToEye, 2*initialDistance) # ideally no far plane, rarg...
    glMatrixMode(GL_MODELVIEW)
    glLoadIdentity()
    glColor3f(.8, 0., 0.)

    ### lighting ###
    glShadeModel (GL_FLAT)
    glEnable(GL_COLOR_MATERIAL)
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, [1.0, 1.0, 1.0, 1.0])
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, [1.0, 1.0, 1.0, 1.0])
    # glMaterialfv(GL_FRONT, GL_SPECULAR, [1.0, 1.0, 1.0, 1.0]) # no specular
    # glMaterialfv(GL_FRONT, GL_SHININESS, [20.0])
    glLightfv(GL_LIGHT0, GL_POSITION, [0, 0.0, initialDistance, 0.0]) # at origin

    glEnable(GL_LIGHTING)
    glEnable(GL_LIGHT0)
    glEnable(GL_DEPTH_TEST)

def square(spin):
    glClearColor(1., 1., 1., 1.)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushMatrix();
    glRotatef(spin, 0.0, 0.0, 1.0);
    glRectf(-25.0, -25.0, 25.0, 25.0);
    glPopMatrix();

# assumes min equivalent to max in a modulus scheme, int args
# forces into range [min, max)
def fix(coord, min, max):
  return min + (coord-min)%(max-min)

def normalize(pt):
  mag = math.sqrt(pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2])
  return [pt[i]/mag for i in range(0, 3)]

# finds cross product of vec<p0, p1> and vec<p0, p2>
def findNormal(p0, p1, p2):
  u = [p1[i]-p0[i] for i in range(0, 3)]
  v = [p2[i]-p0[i] for i in range(0, 3)]
  ret = [0,0,0]
  ret[0] = u[2]*v[1] - u[1]*v[2]
  ret[1] = u[0]*v[2] - u[2]*v[0]
  ret[2] = u[1]*v[0] - u[0]*v[1]
  return normalize(ret)

# each point is a tuple (x,y,z).  Treat first point as anchor
# and create a fan of triangles.  style is int suggesting color
# for normal calculation, we assume all the points are co-planar
def makeConvexPolygon(style, points):
  assert(len(points) >= 3)
  # normal = findNormal(points[0], points[1], points[2])
  # apply(glNormal3f, normal)
  if(style == 0):
    glColor3f(1.0, 0., 0.)
  else:
    glColor3f(0., 0., 1.0)
  glBegin(GL_TRIANGLE_STRIP)
  for p in points:
    # hack - substitute point itself as approximation for its normal
    apply(glNormal3f, normalize(p))
    glVertex3f(p[0], p[1], p[2])
  glEnd()

# lng, lat are in degrees, convert to radians
def convertFromSphereSurfaceToAbsolute(point):
  global radius
  (lng, lat, alt) = point
  r = radius + alt
  y = r * math.sin(lat * 2.0 * math.pi / 360.0)
  minor_radius = math.sqrt(r*r - y*y)
  x = minor_radius * math.cos(lng * 2.0 * math.pi / 360.0)
  z = minor_radius * math.sin(lng * 2.0 * math.pi / 360.0)
  #z = math.sqrt(minor_radius*minor_radius - x*x);
  return (x,y,z)

# each point is a (lng, lat, height from base) on the sphere
# convert to (x,y,z), pass to makeConvexPolygon().  Style is int
# suggesting color.
def makeSphereSurfaceConvexPolygon(style, points):
  new_points = [convertFromSphereSurfaceToAbsolute(p) for p in points]
  makeConvexPolygon(style, new_points)

# on the line of latitude indicated by latLine, span from
# lngStart to lngEnd, by going east or west depending on
# goEast.  The span stretches North/South by width degrees,
# and undershoots and overshoots by width/2 degrees (as though
# a square marker were tracing the span).
def drawLatCorridor(lngStart, lngEnd, latLine, width, goEast=True, granularity=100):

    if(lngStart > lngEnd):
        # just draw the other way.
        drawLatCorridor(lngEnd, lngStart, latLine, width, not goEast, granularity)
        return
    if(not goEast):
        # because of normals and 'crossover', really need to draw three segments.
        drawLatCorridor(-180.0, lngStart, latLine, width, True, granularity)
        drawLatCorridor(-180.0, 180.0, -90 - latLine, width, True, granularity)
        drawLatCorridor(lngEnd, 180.0, latLine, width, True, granularity)
        return

    width /= 2.0 # twice as many degrees latitude as longitude, do this to keep run width same.

    # add under and overshoot, find north and south bounds.
    lngStart -= width/1.0 # not 2 because lng degrees are 1/2 lat degrees.
    lngEnd += width/1.0
    distance = lngEnd-lngStart
    northLat = latLine + width/2.0
    southLat = latLine - width/2.0

    # now, go east from lngStart to lngEnd
    points = []
    for t in xrange(0, granularity+1):
        lng = lngStart + float(t)/granularity * distance
        points.append( (lng, northLat, 0.0) )
        points.append( (lng, southLat, 0.0) )

    # draw
    makeSphereSurfaceConvexPolygon(0, points)

# This is the same as drawLatCorridor, but roles of latitude and
# longitude are reversed.  There's probably some way to generalize.
def drawLngCorridor(latStart, latEnd, lngLine, width, goNorth=True, granularity=100):
    if(latStart > latEnd):
        # just draw the other way.
        drawLngCorridor(latEnd, latStart, lngLine, width, not goNorth, granularity)
        return
    if(not goNorth):
        # because of normals and 'crossover', really need to draw three segments.
        drawLngCorridor(-90.0, latStart, lngLine, width, True, granularity)
        drawLngCorridor(-90.0, 90.0, -180 - lngLine, width, True, granularity)
        drawLngCorridor(latEnd, 90.0, lngLine, width, True, granularity)
        return

    # add under and overshoot, find north and south bounds.
    latStart -= width/4.0 # 4 instead of 2, because lat degrees = 2x lng degrees.
    latEnd += width/4.0
    distance = latEnd-latStart
    westLng = lngLine - width/2.0
    eastLng = lngLine + width/2.0

    # now, go east from latStart to latEnd
    points = []
    for t in xrange(0, granularity+1):
        lat = latStart + float(t)/granularity * distance
        points.append( (westLng, lat, 0.0) )
        points.append( (eastLng, lat, 0.0) )

    # draw
    makeSphereSurfaceConvexPolygon(0, points)

def drawRoutine():
    width = 30.0 - 4.0 # 30.0 max width, 4 degrees of buffer between runs
    drawLngCorridor( -90.0,  60.0,   0.0, width, True)
    drawLatCorridor(   0.0,-150.0,  60.0, width, False)
    drawLngCorridor(  60.0,  30.0,-150.0, width, False)
    drawLatCorridor(-150.0, -30.0,  30.0, width, True)
    drawLngCorridor(  30.0,   0.0, -30.0, width, False)
    drawLatCorridor( -30.0, -60.0,   0.0, width, False)
    drawLngCorridor(   0.0, -30.0, -60.0, width, False)
    drawLatCorridor( -60.0, -30.0, -30.0, width, True)
    drawLngCorridor( -30.0, -60.0, -30.0, width, False)
    drawLatCorridor( -30.0,-150.0, -60.0, width, False)
    drawLngCorridor( -60.0, -30.0,-150.0, width, True)
    drawLatCorridor(-150.0,-120.0, -30.0, width, True)
    drawLngCorridor( -30.0,   0.0,-120.0, width, True)
    drawLatCorridor(-120.0, 120.0,   0.0, width, False)
    drawLngCorridor(   0.0,  30.0, 120.0, width, True)
    drawLatCorridor( 120.0, 150.0,  30.0, width, True)
    drawLngCorridor(  30.0,  60.0, 150.0, width, True)
    drawLatCorridor( 150.0,  60.0,  60.0, width, False)
    drawLngCorridor(  60.0,  30.0,  60.0, width, False)
    drawLatCorridor(  60.0,  90.0,  30.0, width, True)
    drawLngCorridor(  30.0,   0.0,  90.0, width, False)
    drawLatCorridor(  90.0,  60.0,   0.0, width, False)
    drawLngCorridor(   0.0, -30.0,  60.0, width, False)
    drawLatCorridor(  60.0, 150.0, -30.0, width, True)
    drawLngCorridor( -30.0, -60.0, 150.0, width, False)
    drawLatCorridor( 150.0,  30.0, -60.0, width, False)
    drawLngCorridor( -60.0,  90.0,  30.0, width, True)

    # Add polar caps.
    drawLatCorridor(   0.0, 360.0,  80.0, width, True)
    drawLatCorridor(   0.0, 360.0, -80.0, width, True)
    return

def drawCached():
    global displayListHandle
    glCallList(displayListHandle)

def drawInitial():
    global displayListHandle
    print "drawInitial"
    displayListHandle = glGenLists(1)
    print "Gen",displayListHandle
    glNewList(displayListHandle, GL_COMPILE_AND_EXECUTE)
    drawRoutine()
    glEndList()









xrot = 0
yrot = 0
zrot = 0

def keyboardCallback(char, x, y):
    global xrot
    global yrot
    global zrot
    global logDistanceE2
    print char, xrot, yrot, zrot, logDistanceE2
    if(char == 'q'):
        xrot += 1
    if(char == 'z'):
        xrot -= 1
    if(char == 'w'):
        yrot += 1
    if(char == 'x'):
        yrot -= 1
    if(char == 'e'):
        zrot += 1
    if(char == 'c'):
        zrot -= 1
    if(char == 'r'):
        logDistanceE2 += 1
    if(char == 'v'):
        logDistanceE2 -= 1

d = 0
def displayCallback():
    global xrot
    global yrot
    global zrot
    global logDistanceE2
    global d
    global arcBall
    if(d == 0):
        init2()

    glClearColor(0., 0., 0., 0.)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glPushMatrix();
    glMultMatrixf(arcBall.getTransform())
    glTranslatef(0.0,0.0, radius - math.e**(logDistanceE2/100.0));
    glRotatef(xrot, 1.0, 0.0, 0.0)
    glRotatef(yrot, 0.0, 1.0, 0.0)
    glRotatef(zrot, 0.0, 0.0, 1.0)

    if(d == 0):
      drawInitial()
    else:
      drawCached()

    glPopMatrix();
    glFlush();
    glutSwapBuffers()
    glutPostRedisplay()
    d+=1

def init():
    glutInit(name)
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
    glutInitWindowSize(height, width)
    glutCreateWindow(name)
    glutKeyboardFunc(keyboardCallback)
    glutDisplayFunc(displayCallback)
    glutMainLoop()

init()
