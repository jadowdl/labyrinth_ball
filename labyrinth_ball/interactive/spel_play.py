#!/usr/bin/env python2

#########################################
# (c) 2018 James Dowdell, TSM Games
#########################################

"""
It's been almost a decade since I wrote this code.  If I remember right, the
idea was that it is a playpen to explore "spels", or "spherical pixels."
A virtual sphere is grided by lines of latitude and longitude, and each
square cell bordered by two lines of latitude and two lines of longitude
is extruded up a random height, thus causing the sphere to appear hairy.
This ball is then rendered in OpenGL in a window, and manipulated by
keyboard.  It was created mostly to explore the various pyopengl functions.
"""

from OpenGL.GLUT import *
from OpenGL.GLU import *
from OpenGL.GL import *
import random
import math

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

# we want the initial distance to be such that, viewed as a circle on the screen initially
# the circle is centered on the window with radius z, and the window is 4z X 4z.
# Achieved via similar triangles: .5virtsWidth/virtsToEye == 2radius/initialDistance
initialDistance =  2 * radius*virtsToEye / virtsWidth

# distance from surface of sphere, can never get past surface of sphere
logDistanceE2 = int(math.log(initialDistance-radius)*100)


def init2():
  ### generic ###
  glViewport(0,0,400,400)
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


# finds cross product of vec<p0, p1> and vec<p0, p2>
def findNormal(p0, p1, p2):
  u = [p1[i]-p0[i] for i in range(0, 3)]
  v = [p2[i]-p0[i] for i in range(0, 3)]
  ret = [0,0,0]
  ret[0] = u[1]*v[2] - u[2]*v[1]
  ret[1] = u[2]*v[0] - u[0]*v[2]
  ret[2] = u[0]*v[1] - u[1]*v[0]
  mag = math.sqrt(ret[0]*ret[0] + ret[1]*ret[1] + ret[2]*ret[2])
  return [-ret[i]/mag for i in range(0, 3)]


# each point is a tuple (x,y,z).  Treat first point as anchor
# and create a fan of triangles.  style is int suggesting color
# for normal calculation, we assume all the points are co-planar
def makeConvexPolygon(style, points):
  assert(len(points) >= 3)
  normal = findNormal(points[0], points[1], points[2])
  apply(glNormal3f, normal)
  if(style == 0):
    glColor3f(1.0, 0., 0.)
  else:
    glColor3f(0., 0., 1.0)
  glBegin(GL_TRIANGLE_FAN)
  for p in points:
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


# ll = lower left
# ur = upper right
# edgeCase when lllng == urlng, need to define vertices in a different way.
def printConnector(name, style, lllng, lllat, llalt, urlng, urlat, uralt, fixEdge=False): 
  makeSphereSurfaceConvexPolygon(style,
          [(lllng, lllat, llalt),
           (urlng, lllat, llalt if fixEdge else uralt),
           (urlng, urlat, uralt),
           (lllng, urlat, uralt if fixEdge else llalt)])


heights_map = {}
ids_map = {}


def makeSphere():
  global height_maps
  global ids_map
  id = 0
  # a given "point" on the sphere is really a spherical rect
  # (think pixel => "spel").  The coordinate of the spel is the
  # lower left coordinate.
  # We assign random heights to these points.
  for lnge2 in range (-180*prec, 180*prec, 1):
    for late2 in range (-90*prec, 90*prec, 1):
      heights_map[(lnge2, late2)] = radius*.1 *  random.random()
      ids_map[(lnge2, late2)] = id
      id+=1


def drawSphere_():
  global height_maps
  global ids_map
  # Now we need to print all the spels.  Except, we also need to
  # manually extrude them (make them seem like buildings sitting
  # at sea level, whose ceiling is our spel).  an "ExtrusionConnector"
  # is the polygon joining two spels at different heights; rule is that
  # the lower is responsible for adding the connector to the higher.
  # connector should take the style of the higher spel).
  for lnge2 in range (-180*prec,  180*prec, 1):
    for late2 in range (-85*prec, 85*prec, 1):
      id = ids_map[(lnge2, late2)]
      altitude = heights_map[(lnge2, late2)]
      westLngE2 = fix((lnge2-1), -180*prec, 180*prec) # lng of west spel
      eastLngE2 = fix((lnge2+1), -180*prec, 180*prec) # lng of east spel
      northLatE2 = fix((late2+1), -90*prec, 90*prec) # lat of north spel
      southLatE2 = fix((late2-1), -90*prec, 90*prec) # lat of north spel
      westLng = westLngE2*precinv
      eastLng = eastLngE2*precinv
      northLat = northLatE2*precinv
      southLat = southLatE2*precinv
      curLat = late2*precinv
      curLng = lnge2*precinv

      # not strictly a connector, but handy; current spel defined as connector
      # between (curLng, curLat) and (eastLng, northLat), all at same altitude.
      printConnector(id, id%2, 
                     curLng, curLat, altitude, 
                     eastLng, northLat, altitude)
      #printPolygon2(id, id%2,
      #       [curLng,    curLat,   altitude,
      #       curLng,   northLat,   altitude,
      #       eastLng,  northLat,   altitude,
      #       eastLng,    curLat,   altitude,
      #       curLng,     curLat,   altitude])

      # test connect to West
      if( altitude > heights_map[ westLngE2, late2]):
        printConnector("%s-westConnector"%id, id%2,
            curLng, curLat, altitude,
            curLng, northLat, heights_map[ westLngE2, late2])
      # test connect to East
      if( altitude > heights_map[ eastLngE2, late2]):
        printConnector("%s-eastConnector"%id, id%2,
            eastLng, curLat, altitude,
            eastLng, northLat, heights_map[ eastLngE2, late2])
      # test connect to North 
      if( altitude > heights_map[ lnge2, northLatE2]):
        printConnector("%s-northConnector"%id, id%2,
            curLng, northLat, altitude,
            eastLng, northLat, heights_map[ lnge2, northLatE2], True)
      # test connect to South
      if( altitude > heights_map[ lnge2, southLatE2]):
        printConnector("%s-southConnector"%id, id%2,
            curLng, curLat, altitude,
            eastLng, curLat, heights_map[ lnge2, southLatE2], True)


def drawCachedSphere():
  global displayListHandle
  glCallList(displayListHandle)


def drawSphere():
  global displayListHandle
  print "drawSphere"
  displayListHandle = glGenLists(1)
  print "Gen",displayListHandle
  glNewList(displayListHandle, GL_COMPILE_AND_EXECUTE)
  drawSphere_()
  glEndList()


################################################################################

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
  if(d == 0):
    init2()
    makeSphere()
  #print d 

  glClearColor(0., 0., 0., 0.)
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glPushMatrix();
  #glTranslatef(0.0,0.0,-initialDistance + logDistanceE2*10000000.0);
  glTranslatef(0.0,0.0, radius - math.e**(logDistanceE2/100.0));
  #glTranslatef(0.0,0.0, -initialDistance)
  glRotatef(xrot, 1.0, 0.0, 0.0);
  glRotatef(yrot, 0.0, 1.0, 0.0);
  glRotatef(zrot, 0.0, 0.0, 1.0);
  
  if(d == 0):
    drawSphere()
  else:
    drawCachedSphere()
  # glColor3f(0.,1.,0.)
  # glutSolidSphere(radius, 10000, 10000);
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


if __name__ == "__main__":
  init()
