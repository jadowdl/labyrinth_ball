#########################################
# (c) 2011 Evan Mallory, TSM Games
#

from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *

from PIL import Image
import sys
import os
from labyrinth_ball.common.spherical_utils import *
from labyrinth_ball.common.gl_utils import *
from math import pi, cos, sin, sqrt, pow
import numpy
from labyrinth_ball.common.arc_ball import *
import labyrinth_ball.res as res
from scipy.interpolate import interp1d
import random

NO_FROST_POINT =  .93
SO_FROST_POINT = -.93
def setFrostPoint(no, so):
	assert (no >= 0.0 and no <= 1.0)
	assert (so >= -1.0 and so <= 0.0)
	global NO_FROST_POINT
	NO_FROST_POINT = no
	global SO_FROST_POINT
	SO_FROST_POINT = so
	print "set frost points: ", NO_FROST_POINT, SO_FROST_POINT

class xyzInterpolater:
	def __init__(self, maxT, xInterp, yInterp, zInterp):
		self.maxT = maxT
		self.xInterp = xInterp
		self.yInterp = yInterp
		self.zInterp = zInterp

	def value(self, t):
		if t < 0 or t > self.maxT:
			# This is bad when this happens, but it could just be due
			# to a mismatch between the number of green and yellow points,
			# so we will just quietly fail and hope that someone notices.
			return [0, 0, 0]
		else:
			return [self.xInterp(t), self.yInterp(t), self.zInterp(t)]

	def getMaxTValue(self): return self.maxT


class Point:
	KIND_GREEN = 0
	KIND_YELLOW = 1
	KIND_ENDPOINT = 2

	def __init__(self, point, kind, image_point):
		self.point = point
		self.kind = kind
		self.image_point = image_point
		self.neighbors = []
		self.visited = False

	def __repr__(self):
		return "p[%d, %d]" % (self.image_point[0], self.image_point[1])

gPoints = None
gSelectedPoint = None

GREEN = (0, 255, 0)
YELLOW = (255, 255, 0)
PINK = (0xcc, 0, 0xff)

arcBall = None

showPoints = True
showLines = True
fillPath = True
showTrianglePath = False
useInterpolatedPoints = False
numInterps = 1
showBridge = False
bridgeDepth = .6
showImageSphere = False

interpolatedPoints = {}
bridgePoints = {}
bridgeNormals = {}

NORMAL_POINT_SIZE = 4
SELECTED_POINT_SIZE = 8


def getKindForColor(c):
	table = {GREEN : Point.KIND_GREEN,
			 YELLOW : Point.KIND_YELLOW,
			 PINK : Point.KIND_ENDPOINT}
	return table[c]

def getColorForKind(k):
	table = {Point.KIND_GREEN : GREEN,
			 Point.KIND_YELLOW : YELLOW,
			 Point.KIND_ENDPOINT : PINK}
	return table[k]


def drawPoints(kind, points=None):
	global gPoints
	if points is None:
		points = gPoints

	glColor3f(*getColorForKind(kind))
	glBegin(GL_POINTS)
	for p in points:
		if p.kind == kind:
			glVertex3fv(p.point)
	glEnd()


def drawLines(startPoint, kind):
	global interpolatedPoints, useInterpolatedPoints

	if useInterpolatedPoints:
		glEnable(GL_VERTEX_ARRAY)
		glVertexPointer(3, GL_FLOAT, 0, interpolatedPoints[kind].tostring())
		glDrawArrays(GL_LINE_STRIP, 0, len(interpolatedPoints[kind]) / 3)
		glDisable(GL_VERTEX_ARRAY)
	else:
		glBegin(GL_LINE_STRIP)
		while True:
			glVertex3f(*startPoint.point)

			found = False
			for other in startPoint.neighbors:
				if other.kind == kind:
					found = True
					startPoint = other
			if not found:
				break
		glEnd()


def drawPath(curPoint, color, fill=True):
	global interpolatedPoints, useInterpolatedPoints

	glColor4fv(color)
	if not useInterpolatedPoints:
		# We don't have normals if we're not using interpolated points,
		# so don't enable lighting.
		glDisable(GL_LIGHTING)

	if fill:
		glBegin(GL_TRIANGLE_STRIP)
		if not useInterpolatedPoints:
			glVertex3f(*curPoint.point)
	else:
		glBegin(GL_LINES)

	if useInterpolatedPoints:
		for i in range(0, len(interpolatedPoints[Point.KIND_GREEN]) / 3):
			glVertex3fv(interpolatedPoints[Point.KIND_GREEN][i*3:i*3+3])
			glNormal3fv(interpolatedPoints[Point.KIND_GREEN][i*3:i*3+3])
			glVertex3fv(interpolatedPoints[Point.KIND_YELLOW][i*3:i*3+3])
			glNormal3fv(interpolatedPoints[Point.KIND_YELLOW][i*3:i*3+3])
	else:
		one_way = curPoint.neighbors[0]
		the_other_way = curPoint.neighbors[1]
		while True:
			glVertex3fv(one_way.point)
			glVertex3fv(the_other_way.point)

			if not one_way.neighbors or not the_other_way.neighbors:
				break

			one_way = one_way.neighbors[0]
			the_other_way = the_other_way.neighbors[0]
	glEnd()

	if not useInterpolatedPoints:
		glEnable(GL_LIGHTING)

def drawBridge():
	global interpolatedPoints, bridgePoints

	for kind in [Point.KIND_GREEN, Point.KIND_YELLOW]:
		glBegin(GL_TRIANGLE_STRIP)
		for i in range(0, len(interpolatedPoints[kind]) / 3):
			# The upper and lower points on the bridge have the same normal.
			glNormal3fv(bridgeNormals[kind][i*3:i*3+3])
			glVertex3fv(interpolatedPoints[kind][i*3:i*3+3])
			glVertex3fv(bridgePoints[kind][i*3:i*3+3])
		glEnd()


def display():
	global gPoints, arcBall, textureId, sphereModel

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	glLoadIdentity()

	glMultMatrixf(arcBall.getTransform())

	if showPoints:
		glPointSize(NORMAL_POINT_SIZE)
		drawPoints(Point.KIND_GREEN)
		drawPoints(Point.KIND_YELLOW)
		drawPoints(Point.KIND_ENDPOINT)

		if gSelectedPoint:
			glPointSize(SELECTED_POINT_SIZE)
			drawPoints(gSelectedPoint.kind, points=[gSelectedPoint])

	# Draw connecting lines
	curPoint = None
	for p in gPoints:
		if p.kind == Point.KIND_ENDPOINT:
			curPoint = p
			break

	if showLines and curPoint != None:
		glDisable(GL_LIGHTING)
		glColor3f(1, 1, 1)
		drawLines(curPoint, Point.KIND_YELLOW)
		drawLines(curPoint, Point.KIND_GREEN)
		glEnable(GL_LIGHTING)

	if showTrianglePath and curPoint != None:
		drawPath(curPoint, [1, 1, 1, 1], fill=False)
	if fillPath and curPoint != None:
		drawPath(curPoint, [.2, .6, .2, 1], fill=True)

	if showBridge:
		glColor3f(.3, .3, .3)
		drawBridge()

	if showImageSphere:
		glDisable(GL_LIGHTING)
		glColor3f(.6, .6, .6)
		# Make the sphere just a bit smaller than the actual size so that
		# we can see the model that we are editing better.
		glScalef(.97, .97, .97)
		drawModel(sphereModel, texture=textureId)
		glEnable(GL_LIGHTING)

	glutSwapBuffers()


def onKey(key, x, y):
	global fillPath, showTrianglePath, useInterpolatedPoints, numInterps
	global showBridge, bridgeDepth, showPoints, showLines, showImageSphere

	if key == 'f':
		fillPath = not fillPath
		glutPostRedisplay()
	elif key == 't':
		showTrianglePath = not showTrianglePath
		glutPostRedisplay()
	elif key == 'i':
		useInterpolatedPoints = not useInterpolatedPoints
		glutPostRedisplay()
	elif key == 'p':
		showPoints = not showPoints
		glutPostRedisplay()
	elif key == 'l':
		showLines = not showLines
		glutPostRedisplay()
	elif key == 'b':
		showBridge = not showBridge
		glutPostRedisplay()
	elif key == 's':
		showImageSphere = not showImageSphere
		glutPostRedisplay()
	elif key == '1': # wireframe
		showPoints = True
		showLines = True
		fillPath = False
		showTrianglePath = False
		useInterpolatedPoints = False
		showBridge = False
		showImageSphere = True
		glutPostRedisplay()
	elif key == '2': # solid
		showPoints = False
		showLines = False
		showTrianglePath = False
		fillPath = True
		useInterpolatedPoints = True
		showBridge = True
		showImageSphere = False
		glutPostRedisplay()
	elif key == 'e':
		exportModel()
	elif key == 'q':
		writeMapImage()
		sys.exit(0)
	elif key == 'w':
		writeMapImage()
	elif key == '=': # +
		numInterps += 1
		interpolatePoints()
		glutPostRedisplay()
	elif key == '-':
		if numInterps > 1:
			numInterps -= 1
		interpolatePoints()
		glutPostRedisplay()
	elif key == ',': # <
		bridgeDepth -= .1
		buildBridge()
		glutPostRedisplay()
	elif key == '.': # >
		bridgeDepth += .1
		buildBridge()
		glutPostRedisplay()
	elif key == 'r':
		loadMapImage()
		loadSinusoidalMapAsTexture(sys.argv[1])
		glutPostRedisplay()


def onSpecialKey(key, x, y):
	if gSelectedPoint is None:
		return

	(_, lat, lng) = xyzToSphericalCoords(gSelectedPoint.point)

	if key == GLUT_KEY_RIGHT:
		lng += pi / 180
	elif key == GLUT_KEY_LEFT:
		lng -= pi / 180
	elif key == GLUT_KEY_UP:
		lat -= pi / 180
	elif key == GLUT_KEY_DOWN:
		lat += pi / 180

	if lng < -pi:
		lng += 2 * pi
	if lng > pi:
		lng -= 2 * pi
	if lat < 0:
		lat = 0
	if lat > pi:
		lat = pi

	gSelectedPoint.point = sphericalCoordsToXyz(1, lat, lng)
	interpolatePoints()
	glutPostRedisplay()


def writeVertex(file, point, normal):
	file.write("%f %f %f %f %f %f\n" % (point[0], point[1], point[2],
		normal[0], normal[1], normal[2]))


def writeVertexInTriangleStrip(file, vertex_indices, point, normal, vertex_i):
	writeVertex(file, point, normal)

	vertex_indices.append(vertex_i)


def averagePoints(p1, p2, p1_weight):
	average = [0, 0, 0]
	for i in range(0, 3):
		average[i] = p1_weight * p1[i] + (1 - p1_weight) * p2[i]
	return normalizePoint(average)


def exportModel():
	global NO_FROST_POINT, SO_FROST_POINT
	out_filename = os.path.splitext(sys.argv[1])[0] + ".txt"
	file = open(out_filename, 'w')

	# Write texture coordinates
	file.write('pre_computed\n')
	# texture coordinates for roof
	pointsPerLength = 20
	for i in range(0, len(interpolatedPoints[Point.KIND_GREEN]) / 3):
		# Normally we want the texture coordinate for the roof perpendicular
		# to run from .2 to .8, but if we are close enough to the
		# pole at (0, 0, +/-1), limit the range to .3 to .7.  Clients
		# should update the frost point for their pole by calling
		# setFrostPoint().  This is  so that the white "frosting"
		# at the edge of the path isn't shown where we don't want it to be.
		if (interpolatedPoints[Point.KIND_GREEN][i*3 + 2] > NO_FROST_POINT or
		   interpolatedPoints[Point.KIND_GREEN][i*3 + 2] < SO_FROST_POINT):
			min_val = .3
		else:
			min_val = .2
		if (interpolatedPoints[Point.KIND_YELLOW][i*3 + 2] > NO_FROST_POINT or
		   interpolatedPoints[Point.KIND_YELLOW][i*3 + 2] < SO_FROST_POINT):
			max_val = .7
		else:
			max_val = .8

		for j in [1, 0]:
			distToEdge = (float(j) * (max_val - min_val)) + min_val
			file.write("%f %f\n" % (distToEdge,
				float(i) / numInterps / pointsPerLength))

	# texture coordinates for bridge
	for kind in (Point.KIND_GREEN, Point.KIND_YELLOW):
		for i in range(0, len(interpolatedPoints[kind]) / 3):
			# two points per interpolated point. The upper point comes first.

			v = float(i) / numInterps / pointsPerLength

			if kind == Point.KIND_GREEN:
				file.write("-.2 %f\n" % v)
				file.write(".2 %f\n" % v)
			else:
				file.write(".8 %f\n" % v)
				file.write("1.2 %f\n" % v)

	vertex_indices = []
	vertex_i = 0

	# Write roof
	file.write('\n')
	for i in range(0, len(interpolatedPoints[Point.KIND_GREEN]) / 3):
		p1 = interpolatedPoints[Point.KIND_GREEN][i*3:i*3+3]
		p2 = interpolatedPoints[Point.KIND_YELLOW][i*3:i*3+3]

		for j in [1, 0]:
			p3 = averagePoints(p1, p2, j)

			writeVertexInTriangleStrip(file,
				vertex_indices,
				p3,
				p3,
				vertex_i)

			vertex_i += 1

	# Write bridge
	for kind in (Point.KIND_GREEN, Point.KIND_YELLOW):
		for i in range(0, len(interpolatedPoints[kind]) / 3):
			# TODO: Don't connect the first bridge points with the
			# last of the path points from above or the previous color's
			# bridge points

			# For one side of the bridge, we write the upper vertex first,
			# and for the other, we write the lower vertex first so that
			# all the triangles are listed in clockwise order.
			if kind == Point.KIND_YELLOW:
				writeVertexInTriangleStrip(file,
					vertex_indices,
					interpolatedPoints[kind][i*3:i*3+3],
					bridgeNormals[kind][i*3:i*3+3],
					vertex_i)
				vertex_i += 1

			writeVertexInTriangleStrip(file,
				vertex_indices,
				bridgePoints[kind][i*3:i*3+3],
				bridgeNormals[kind][i*3:i*3+3],
				vertex_i)
			vertex_i += 1

			if kind == Point.KIND_GREEN:
				writeVertexInTriangleStrip(file,
					vertex_indices,
					interpolatedPoints[kind][i*3:i*3+3],
					bridgeNormals[kind][i*3:i*3+3],
					vertex_i)
				vertex_i += 1

	file.write("\n")

	file.write("TRIANGLE_STRIP\n\n")

	for i in vertex_indices:
		file.write("%d\n" % i)

	file.write("\n")

	# Write path points
	curPoint = None
	for p in gPoints:
		if p.kind == Point.KIND_ENDPOINT:
			curPoint = p
			break

	one_way = curPoint.neighbors[0]
	the_other_way = curPoint.neighbors[1]
	num_path_points = 0
	path_points = []
	north_pole_sucked = False
	south_pole_sucked = False
	while True:
		num_path_points += 1
		average = averagePoints(one_way.point, the_other_way.point, .5)
		# Suck path points above/below the frost point to the pole, and don't
		# duplicate.  Useful for both manual models (hard to deal with sinuisoidal
		# near the pole) and autogenerated (which spirals around at the pole).
		if (average[2] > NO_FROST_POINT):
			if (not north_pole_sucked):
				north_pole_sucked = True
				path_points.append([0,0,1])
		elif (average[2] < SO_FROST_POINT):
			if (not south_pole_sucked):
				south_pole_sucked = True
				path_points.append([0,0,-1])
		else:
			path_points.append(average)

		if not one_way.neighbors or not the_other_way.neighbors:
			break

		one_way = one_way.neighbors[0]
		the_other_way = the_other_way.neighbors[0]
	# Hack - because of our cubic interpolation, and the sucking action above,
	# if we leave the poles unattended we might get weird kinks in the path
	# going to/from the pole.  Mitigate via manual linear interpolation.
	pole_fix_count = 5
	if(path_points[0][2] == 1.0 or path_points[0][2] == -1.0):
		a = path_points[0]
		b = path_points[1]
		for i in range(1, pole_fix_count):
		    path_points.insert(1, averagePoints(a, b, float(i)/pole_fix_count))
        # make the pole point just off center to avoid mathematical singularity
        path_points[0] = averagePoints(path_points[0], path_points[1], .99)
	if(path_points[-1][2] == 1.0 or path_points[-1][2] == -1.0):
		a = path_points[-1]
		b = path_points[-2]
		for i in range(1, pole_fix_count):
		    path_points.insert(-1, averagePoints(a, b, float(i)/pole_fix_count))
        # make the pole point just off center to avoid mathematical singularity
        path_points[-1] = averagePoints(path_points[-1], path_points[-2], .99)

	# Actually write the path.
	for p in path_points:
		file.write("%f %f %f 0\n" % (p[0], p[1], p[2]))

	file.close()
	print "Exported to %s" % out_filename
	print "\t%d vertices, %d vertex indices, %d path points" % (
		vertex_i, len(vertex_indices), num_path_points)


def onResize(width, height):
	glViewport(0, 0, width, height)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	ratio = float(width) / height
	small_side = 1.05
	large_side = small_side * max(width, height) / min(width, height)
	if width > height:
		glOrtho(-large_side, large_side, -small_side, small_side, -2, .1)
	else:
		glOrtho(-small_side, small_side, -large_side, large_side, -2, .1)
	glMatrixMode(GL_MODELVIEW)
	glLoadIdentity()
	arcBall.setBounds(width, height)


def loadSinusoidalMapAsTexture(filename):
	global textureId

	textureId = glGenTextures(1)

	glBindTexture(GL_TEXTURE_2D, textureId)
	interpMode = GL_LINEAR
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, interpMode)
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, interpMode)

	image = Image.open(filename)
	(width, height) = image.size
	image = image.convert("RGB")

	image_data = image.tobytes("raw", "RGB", 0, -1)
	glTexImage2D(GL_TEXTURE_2D, 0, 3, width, height, 0, GL_RGB,
		GL_UNSIGNED_BYTE, image_data)


def onClick(button, button_state, x, y):
	global gSelectedPoint, gPoints

	if button == GLUT_LEFT_BUTTON:
		arcBall.onClick(button, button_state, x, y)
	elif button == GLUT_RIGHT_BUTTON:
		# Select the point nearest the click
		bestDist = 100000
		bestPoint = None

		clickPoint = arcBall.mapToSphere(Point2f(x, y))
		spherePoint = numpy.dot(arcBall.getTransform()[:3,:3], clickPoint)

		for p in gPoints:
			dist = sqrt(pow(p.point[0] - spherePoint[0], 2) +
				pow(p.point[1] - spherePoint[1], 2) +
				pow(p.point[2] - spherePoint[2], 2))
			if dist < bestDist:
				bestPoint = p
				bestDist = dist

		gSelectedPoint = bestPoint
		glutPostRedisplay()


def displayPoints():
	global arcBall, sphereModel

	glutInit(sys.argv)
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH)
	glutInitWindowSize(400, 400)
	arcBall = ArcBall(400, 400)
	glutInitWindowPosition(100, 100)
	glutCreateWindow(sys.argv[0])

	glClearColor(.3, .5, .7, 1)
	glEnable(GL_DEPTH_TEST)
	glShadeModel(GL_SMOOTH)

	glEnable(GL_LIGHTING)
	glEnable(GL_LIGHT0)
	glLightfv(GL_LIGHT0, GL_POSITION, [1, 1, 3, 0])
	glLightfv(GL_LIGHT0, GL_AMBIENT, [0, 0, 0, 1])
	glEnable(GL_COLOR_MATERIAL)

	glutDisplayFunc(display)
	glutMouseFunc(onClick)
	glutMotionFunc(lambda x, y: arcBall.onMouseMove(x, y))
	glutReshapeFunc(onResize)
	glutKeyboardFunc(onKey)
	glutSpecialFunc(onSpecialKey)

	loadSinusoidalMapAsTexture(sys.argv[1])
        resDir = os.path.abspath(res.__path__[0]) 
	sphereModel = loadModel(os.path.join(resDir, "sinusoidal_sphere.3d"))

	glutMainLoop()


def loadPoints(image):
	(width, height) = image.size
	pixels = image.load()
	result = []
	for y in range(0, height):
		for x in range(0, width):
			color = pixels[x, y]
			if (color == GREEN or
				color == YELLOW or
				color == PINK):
				# compute lat/long
				(latitude, longitude) = sinusoidal2dToLatLong(
					x, y, width, height)

				if (latitude <= pi/2 and latitude >= -pi/2 and
					longitude <= pi and longitude >= -pi):
					latitude += pi/2

					point = numpy.array([sin(latitude) * cos(longitude),
										 sin(latitude) * sin(longitude),
										 cos(latitude)], 'f')

					result.append(Point(point, getKindForColor(color), (x,y)))
	return result


def loadMapImage():
	global gPoints

	image = Image.open(sys.argv[1])
	if image.mode != 'RGB':
		image = image.convert("RGB")
	(width, height) = image.size
	assert width == height * 2

	print "Loading points"
	gPoints = loadPoints(image)
	print "Computing neighbors"
	if(len(gPoints) > 0):
		computePointNeighbors(gPoints, image)
		interpolatePoints()

def loadPrecomputedPoints(points):
	global gPoints
	print "Updating gPoints"
	gPoints = points
	if(len(gPoints) > 0):
		interpolatePoints()

def writeMapImage():
	global gPoints

	print "Writing map image..."

	image = Image.open(sys.argv[1])
	if image.mode != 'RGB':
		image = image.convert("RGB")
	(width, height) = image.size
	pixels = image.load()

	out_image = Image.new("RGB", (width, height), (255, 255, 255))
	out_pixels = out_image.load()

	# Write everything but the point locations
	for y in range(0, height):
		for x in range(0, width):
			in_pixel = pixels[x, y]

			(latitude, longitude) = sinusoidal2dToLatLong(
				x, y, width, height)

			# If a pixel in the range of the actual map is a point color
			# pixel, write the background color. Otherwise, write the
			# original pixel.
			if (latitude <= pi/2 and latitude >= -pi/2 and
				longitude <= pi and longitude >= -pi and
				(in_pixel == GREEN or in_pixel == YELLOW or
					in_pixel == PINK)):
				out_pixels[x, y] = (70, 70, 70)
			else:
				out_pixels[x, y] = in_pixel

	# Write the point locations
	for p in gPoints:
		(x, y) = xyzToSinusoidal2d(p.point, width, height)
		out_pixels[x, y] = getColorForKind(p.kind)

	out_image.save("out-" + sys.argv[1])



def dist(p1, p2):
	# spherical distance between two unit vectors
	return acos(numpy.dot(p1, p2))
	# return sqrt(sum([c*c for c in p2-p1]))


# Return True iff you can get from p1 to p2 without
# crossing a black line on the image.
def canSee(image, p1, p2, debug=False):
	dx = p2.image_point[0] - p1.image_point[0]
	dy = p2.image_point[1] - p1.image_point[1]
	num_steps = max(abs(dx), abs(dy))

	(width, height) = image.size
	pixels = image.load()

	dx = p2.point[0] - p1.point[0]
	dy = p2.point[1] - p1.point[1]
	dz = p2.point[2] - p1.point[2]
	if debug: print "canSee from %s to %s" % (p1, p2)

	for step in range(0, num_steps + 1):
		x2 = p1.point[0] + (float(step) / num_steps) * dx
		y2 = p1.point[1] + (float(step) / num_steps) * dy
		z2 = p1.point[2] + (float(step) / num_steps) * dz
		(img_x, img_y) = xyzToSinusoidal2d([x2, y2, z2], width, height)
		if pixels[img_x, img_y] == (0,0,0):
			if debug: print "Found wall at %d, %d (%.2f, %.2f, %.2f)" % (
				img_x, img_y,
				x2, y2, z2)
			return False
	return True


def cmpPointDist(x, y):
	if x[1] < y[1]:
		return -1
	elif x[1] == y[1]:
		return 0
	else:
		return 1


def clearVisited(points):
	for p in points:
		p.visited = False


def computePointNeighborsFromPoint(points, image, curPoint, kind):
	clearVisited(points)

	while True:
		# array of (Point, distance)
		curPoint.visited = True
		distances = []

		for other in points:
			if (curPoint != other and
				not other.visited and
				other.kind == kind):
				distances.append((other, dist(curPoint.point, other.point)))

		# sort by distance
		distances.sort(cmpPointDist)

		foundOther = False
		if len(distances) > 0:
			# Add closest point as a neighbor if we can see it.
			# Defer canSee calculation to here because it is slow.
			for candidate in distances:
				other = candidate[0]
				if canSee(image, curPoint, other):
					curPoint.neighbors.append(other)
					curPoint = other
					foundOther = True
					break
		if not foundOther:
			break

	hasUnconnectedPoints = False
	connectedPointCount = 0
	for p in points:
		if p.kind == kind:
			if p.visited == False:
				if not hasUnconnectedPoints:
					hasUnconnectedPoints = True
					print "\tUnconnected points of kind %d: " % (kind),
				print p,
			else:
				connectedPointCount += 1
	if hasUnconnectedPoints:
		print
	print "\t%d connected points" % connectedPointCount


def computePointNeighbors(points, image):
	curPoint = None
	for p in points:
		if p.kind == Point.KIND_ENDPOINT:
			curPoint = p
			break

	computePointNeighborsFromPoint(points, image, curPoint, Point.KIND_GREEN)
	computePointNeighborsFromPoint(points, image, curPoint, Point.KIND_YELLOW)


def getInterpolater(curPoint, points, kind, useActualInterpolater):
	xvals = []
	yvals = []
	zvals = []
	tvals = []

	numPoints = 0
	keepLooping = True
	while keepLooping:
		tvals.append(numPoints)

		# Snap points to the poles if they are close enough.
		(_, lat, lng) = xyzToSphericalCoords(curPoint.point)
		if lat > (pi * .95):
			lat = pi
			(x, y, z) = sphericalCoordsToXyz(1, lat, lng)
		elif lat < (pi * .05):
			lat = 0
			(x, y, z) = sphericalCoordsToXyz(1, lat, lng)
		else:
			(x, y, z) = curPoint.point

		xvals.append(x)
		yvals.append(y)
		zvals.append(z)
		numPoints += 1

		keepLooping = False
		for p in curPoint.neighbors:
			if p.kind == kind:
				keepLooping = True
				curPoint = p
				break

	interpKind = 'cubic'
	if useActualInterpolater:
		return xyzInterpolater(
			numPoints - 1,
			interp1d(tvals, xvals, kind=interpKind),
			interp1d(tvals, yvals, kind=interpKind),
			interp1d(tvals, zvals, kind=interpKind))
	else:
		return xyzInterpolater(
			numPoints - 1,
			no_interp1d(tvals, xvals, kind=interpKind),
			no_interp1d(tvals, yvals, kind=interpKind),
			no_interp1d(tvals, zvals, kind=interpKind))

class NoInterp1d:
	def __init__(self, vals):
		self.vals = vals

	def __call__(self, t):
		i = int(t)
		return self.vals[i]

def no_interp1d(tvals, xvals, kind):
	return NoInterp1d(xvals)

def computeNormal(p1, p2, p3, kind, debug=False):
	p1 = numpy.array(p1, 'f')
	p2 = numpy.array(p2, 'f')
	p3 = numpy.array(p3, 'f')

	v1 = p1 - p2
	v2 = p3 - p2

	mult = 1
	# TODO: fix this heinous hack
	if kind == Point.KIND_GREEN:
		mult = -1
	return mult * normalizePoint(numpy.cross(v1, v2))


def buildBridge():
	global interpolatedPoints, gPoints, bridgePoints, bridgeDepth, bridgeNormals

	print "Building bridge of depth %.2f" % bridgeDepth

	for kind in [Point.KIND_GREEN, Point.KIND_YELLOW]:
		bridgePoints[kind] = numpy.zeros(len(interpolatedPoints[kind]), 'f')
		bridgeNormals[kind] = numpy.zeros(len(interpolatedPoints[kind]), 'f')
		for i in range(0, len(bridgePoints[kind]) / 3):
			for j in range(0, 3):
				bridgePoints[kind][i*3+j] = (
					interpolatedPoints[kind][i*3+j] * bridgeDepth)

			if i > 0:
				di = -1
			else:
				di = 1
			normal = computeNormal(bridgePoints[kind][i*3:i*3+3],
				interpolatedPoints[kind][i*3:i*3+3],
				interpolatedPoints[kind][(i+di)*3:(i+di)*3+3],
				kind)

			# budge the bridge out a bit in the normal direction
			for j in range(0, 3):
				bridgePoints[kind][i*3+j] += normal[j] * .15

			# Then recompute the normal based on the new location
			normal = computeNormal(bridgePoints[kind][i*3:i*3+3],
				interpolatedPoints[kind][i*3:i*3+3],
				interpolatedPoints[kind][(i+di)*3:(i+di)*3+3],
				kind)

			for j in range(0, 3):
				bridgeNormals[kind][i*3+j] = normal[j]


def interpolatePoints():
	global interpolatedPoints, numInterps, gPoints

	print "Interpolating points with %d divisions." % numInterps

	startPoint = None
	for p in gPoints:
		if p.kind == Point.KIND_ENDPOINT:
			startPoint = p
			break

	for kind in [Point.KIND_GREEN, Point.KIND_YELLOW]:
		interpolater = getInterpolater(startPoint, gPoints, kind, numInterps > 1)
		interpolatedPoints[kind] = numpy.zeros(
			3 * (interpolater.getMaxTValue() * numInterps + 1), 'f')

		print "\t%d vertices" % (len(interpolatedPoints[kind]) / 3)

		for t in range(0, len(interpolatedPoints[kind]) / 3):
			val = interpolater.value(t / float(numInterps))
			val = normalizePoint(val)
			for i in range(0, 3):
				interpolatedPoints[kind][t*3 + i] = val[i]

	# Points starting from 2 before the end of the green should be
	# the same as points starting 21 before the end of the yellow.
	# Total hack to get labyrinth_7.png model just right at the north pole.
	#greenStart = 2 * numInterps
	#yellowStart = 21 * numInterps
	#numPoints = len(interpolatedPoints[Point.KIND_GREEN]) / 3
	#for i in range(numPoints - greenStart, numPoints):
	#	for j in range(0, 3):
	#		interpolatedPoints[Point.KIND_GREEN][i*3 + j] = (
	#			interpolatedPoints[Point.KIND_YELLOW]
	#			[(i + greenStart - yellowStart)*3 + j])

	buildBridge()


def main ():
	assert len(sys.argv) == 2,\
		"Usage: %s <input_image>" % sys.argv[0]

	loadMapImage()

	displayPoints()


if __name__ == "__main__":
	main()

