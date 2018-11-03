#########################################
# (c) 2011 Evan Mallory, TSM Games
#

from utils import *

from OpenGL.GLUT import *
from OpenGL.GL import *
from OpenGL.GLU import *

from array import array
import sys
from time import sleep

from ArcBall import *
import numpy
from scipy.interpolate import interp1d
from math import sin, cos, pi, acos, sqrt

# From "PIL", the Python Imaging Library
import Image

mainModel = None
waypointModel = None

textureIds = None

arcBall = None

showPath = True
showWaypoints = True

PATH_MOVE_INCREMENT = .7
PATH_CLICK_LOOKAHEAD = PATH_MOVE_INCREMENT * 20
PATH_DRAW_INCREMENT = .1

pathDirection = 1 # -1 if we're going backwards
curPathPointT = 0

# interpolation between path point functions
xInterp = None
yInterp = None
zInterp = None

def extrude(point, multiplier):
	return [c * multiplier for c in point]


def pathPointAt(t):
	global xInterp, yInterp, zInterp

	return [xInterp(t), yInterp(t), zInterp(t)]


def display():
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
	glLoadIdentity()
	glEnable(GL_COLOR_MATERIAL)

	glPushMatrix()
	glMultMatrixf(arcBall.getTransform())

	# Draw the model
	glEnable(GL_LIGHTING)
	glColor3f(1, 1, 1)
	drawModel(mainModel, texture=textureIds[0])

	# Disable lighting so we don't have to worry about normals.
	glDisable(GL_LIGHTING)
	if mainModel.path_points:
		denom = 255.0
		# colors = [12/denom, 15/denom, 50/denom, 1]
		colors = [146/denom, 183/denom, 207/denom, 1]
		glColor4fv(colors)

		# Fancy line code
#   	glEnable(GL_LINE_SMOOTH)
#   	glEnable(GL_BLEND)
#   	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA)
#   	glLineWidth(2)

		glBegin(GL_LINE_STRIP)
		t = 0
		while t < maxPathDistance:
			t += PATH_DRAW_INCREMENT
			point = normalizePoint(pathPointAt(t))
			if showPath:
				glVertex3fv(point * 1.01)

		glEnd()

#   	glDisable(GL_BLEND)
#   	glDisable(GL_LINE_SMOOTH)

	# Draw path points from main model, remembering waypoints & goals
	waypoints = []
	goals = []
	for point in mainModel.path_points:
		if point[3] == 1:
			waypoints.append(point[0:3])
		if i == len(mainModel.path_points) - 1:
			goals.append(point[0:3])

	glEnable(GL_LIGHTING)
	glColor3f(1, 1, 1)
	# glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, [.5, .5, .5, 1]);
	if waypointModel and showWaypoints:
		# Draw waypoint model at each waypoint
		for point in waypoints:
			glPushMatrix()
			(r, theta, phi) = xyzToSphericalCoords(point)

			glRotatef(phi * 180/pi, 0, 0, 1)
			glRotatef(theta * 180/pi, 0, 1, 0)
			glTranslate(0, 0, 1)

			if len(textureIds) > 1:
				drawModel(waypointModel, texture=textureIds[1])
			else:
				glColor3f(1, 1, 0)
				drawModel(waypointModel)
				glColor3f(1, 1, 1)

			glPopMatrix()

	# Draw waypoint model at current point
#   if mainModel.path_points:
#   	glTranslate(*pathPointAt(curPathPointT))
#   	glRotatef(90, 1, 0, 0)
#   	drawModel(waypointModel, texture=textureIds[1])

	# Reset emissive to zero.
	glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, [0, 0, 0, 1]);

	glPopMatrix()
	glutSwapBuffers()


def loadTextures(filenames):
	global textureIds

	textureIds = glGenTextures(len(filenames))

	# If we only ask for one texture, it doesn't get put into an array,
	# so we do it ourselves here.
	if len(filenames) == 1:
		textureIds = [textureIds]

	for i in range(0, len(filenames)):
		glBindTexture(GL_TEXTURE_2D, textureIds[i])
		interpMode = GL_NEAREST
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, interpMode)
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, interpMode)

		image = Image.open(filenames[i])
		(width, height) = image.size
		image = image.convert("RGBX")

		image_data = image.tostring("raw", "RGBX", 0, -1)
		glTexImage2D(GL_TEXTURE_2D, 0, 3, width, height, 0, GL_RGBA,
			GL_UNSIGNED_BYTE, image_data)


def onResize(width, height):
	global gScreenWidth, gScreenHeight, gLogicalWidth, gLogicalHeight

	glViewport(0, 0, width, height)
	glMatrixMode(GL_PROJECTION)
	glLoadIdentity()
	ratio = float(width) / height
	small_side = 1.0
	large_side = small_side * max(width, height) / min(width, height)
	if width > height:
		gLogicalWidth = 2 * large_side
		gLogicalHeight = 2 * small_side
		glOrtho(-large_side, large_side, -small_side, small_side, -2, 2)
	else:
		gLogicalWidth = 2 * small_side
		gLogicalHeight = 2 * large_side
		glOrtho(-small_side, small_side, -large_side, large_side, -2, 2)
	glMatrixMode(GL_MODELVIEW)
	glLoadIdentity()
	arcBall.setBounds(width, height)

	gScreenWidth = width
	gScreenHeight = height


def rotateToCurPathPoint():
	global arcBall

	if not mainModel.path_points:
		return

	rot = numpy.identity(3, 'f')
	trans = numpy.identity(4, 'f')

	# Don't calculate as if we are at the absolute end (even if we are),
	# because the rotation calculation below goes wacky.
	point = pathPointAt(min(curPathPointT, maxPathDistance - .01))
	(r, theta, phi) = xyzToSphericalCoords(point)

	# Let OpenGL do that multiplication for us
#   glMatrixMode(GL_MODELVIEW)
#   glLoadIdentity()
#   glRotatef(-theta * 180 / pi, 0, 1, 0)
#   glRotatef(-phi * 180 / pi, 0, 0, 1)
#   trans = glGetDoublev(GL_MODELVIEW_MATRIX)

	# Do the math ourselves
	zrot = numpy.identity(3, 'f')
	# phi = -phi
	zrot[0][0] = cos(phi)
	zrot[0][1] = -sin(phi)
	zrot[1][0] = sin(phi)
	zrot[1][1] = cos(phi)

	theta = -theta
	yrot = numpy.identity(3, 'f')
	yrot[0][0] = cos(theta)
	yrot[0][2] = -sin(theta)
	yrot[2][0] = sin(theta)
	yrot[2][2] = cos(theta)

	rot = numpy.dot(zrot, yrot)

	# Rotate the "North Pole" based on our computed rotation matrix
	rotatedPoint = numpy.dot([0, 0, 1], rot)
	# We want the pole to be position at 90 degrees phi ("up").
	# So we need to rotate the difference between the current
	# angle and the desired one.
	phi = xyzToSphericalCoords(rotatedPoint)[2] - pi/2

	zrot2 = numpy.identity(3, 'f')
	zrot2[0][0] = cos(phi)
	zrot2[0][1] = -sin(phi)
	zrot2[1][0] = sin(phi)
	zrot2[1][1] = cos(phi)

	rot = numpy.dot(rot, zrot2)

	arcBall.setRotation(rot)
	glutPostRedisplay()


def moveTowardsPathPoint(goBackwards):
	global curPathPointT, pathDirection

	direction = pathDirection
	if goBackwards:
		direction = -direction

	curPathPointT += PATH_MOVE_INCREMENT * direction
	curPathPointT = max(0, min(maxPathDistance, curPathPointT))

	rotateToCurPathPoint()


def onKey(key, x, y):
	global showPath, showWaypoints, pathDirection

	if key == 'p':
		showPath = not showPath
		glutPostRedisplay()
	elif key == 'w':
		showWaypoints = not showWaypoints
		glutPostRedisplay()
	elif key == 'j':
		moveTowardsPathPoint(False)
	elif key == 'k':
		moveTowardsPathPoint(True)
	elif key == 'q':
		sys.exit(0)
	elif key == 'r':
		reloadAllFiles()
		glutPostRedisplay()


def onClick(button, button_state, x, y):
	global arcBall, curPathPointT, maxPathDistance
	global gScreenWidth, gScreenHeight, gLogicalWidth, gLogicalHeight

	if button == GLUT_LEFT_BUTTON:
		arcBall.onClick(button, button_state, x, y)
	elif button == GLUT_RIGHT_BUTTON and button_state == GLUT_DOWN:
		# Position marker close to the click
		point = [0, 0]
		point[0] = (x / (gScreenWidth - 1.0)) * gLogicalWidth - (gLogicalWidth / 2.0)
		point[1] = - (y / (gScreenHeight - 1.0)) * gLogicalHeight + (gLogicalHeight / 2.0)

		dist_sq = numpy.dot(point, point)
		if dist_sq < .9: # Don't place based on clicks too near the edge
			click_point = [0, 0, 0]
			click_point[0] = point[0]
			click_point[1] = point[1]
			click_point[2] = sqrt(1 - dist_sq)

			# Rotate click_point based on our current rotation matrix
			transform = arcBall.getTransform()
			rotation = numpy.identity(3, 'f')
			for i in range(0, 3):
				for j in range(0, 3):
					rotation[i][j] = transform[i][j]
			sphere_point = numpy.dot(rotation, click_point)

			# Find the distance from points on the path near the current
			# location to this new point. If the closest one is close enough,
			# set it to be the current point.
			closestT = -1
			closestDistance = 1000000
			for t in drange(max(0, curPathPointT - PATH_CLICK_LOOKAHEAD),
				min(maxPathDistance, curPathPointT + PATH_CLICK_LOOKAHEAD),
				PATH_CLICK_LOOKAHEAD / 50.0):
				p = pathPointAt(t)
				d = dist(p, sphere_point)
				if d < closestDistance:
					closestDistance = d
					closestT = t

			if closestDistance < .1:
				curPathPointT = closestT
				rotateToCurPathPoint()


def dist(p1, p2):
	a1 = normalizePoint(p1)
	a2 = normalizePoint(p2)
	return sqrt(sum([x*x for x in a2-a1]))


def interp1d_wrapper(xs, path_points, point_index):
	ys = [p[point_index] for p in path_points]
	return interp1d(xs, ys, bounds_error = False, fill_value = ys[-1], kind='cubic')

def interpolatePathPoints(model):
	global xInterp, yInterp, zInterp, maxPathDistance

	if not model.path_points:
		return

	distances = numpy.array([0], 'f')

	for i in range(1, len(model.path_points)):
		distances = numpy.append(distances,
			distances[-1] + dist(model.path_points[i-1][:3],
				model.path_points[i][:3]))

	maxPathDistance = distances[-1]

	xInterp = interp1d_wrapper(distances, model.path_points, 0)
	yInterp = interp1d_wrapper(distances, model.path_points, 1)
	zInterp = interp1d_wrapper(distances, model.path_points, 2)


def reloadAllFiles():
	global mainModel, waypointModel

	mainModel = loadModel(sys.argv[1])
	if len(sys.argv) > 3:
		waypointModel = loadModel(sys.argv[3])

	interpolatePathPoints(mainModel)

	if len(sys.argv) > 4:
		loadTextures([sys.argv[2], sys.argv[4]])
	else:
		loadTextures([sys.argv[2]])

	glBindTexture(GL_TEXTURE_2D, textureIds[0])
	if mainModel.clampS:
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE)
	if mainModel.clampT:
		glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE)

def main():
	assert len(sys.argv) >= 3 and len(sys.argv) <= 5, \
		"Usage: %s <model_file> <texture_file> " \
		"[<waypoint_file> [<waypoint texture>]]" % sys.argv[0]

	global arcBall, curPathPointT, mainModel

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
	glLightfv(GL_LIGHT0, GL_POSITION, [1, 1, 3, 1])
	glLightfv(GL_LIGHT0, GL_AMBIENT, [0, 0, 0, 1])
	glLightfv(GL_LIGHT0, GL_DIFFUSE, [.4, .4, .4, 1])

	glEnable(GL_LIGHT1)
	glLightfv(GL_LIGHT1, GL_POSITION, [0, 0, 2, 1])
	glLightfv(GL_LIGHT1, GL_AMBIENT, [0, 0, 0, 1])
	glLightfv(GL_LIGHT1, GL_DIFFUSE, [1, 1, 1, 1])
	glLightf(GL_LIGHT1, GL_SPOT_EXPONENT, 50)
	glLightf(GL_LIGHT1, GL_SPOT_CUTOFF, 30)
	glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, [0, 0, -1])

	glFrontFace(GL_CW)
	glCullFace(GL_BACK)
	glEnable(GL_CULL_FACE)

	# This is needed to see specular highlights on top of textures.
	# I think it is not available on Android, unfortunately.
	# glLightModelfv(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR)

	reloadAllFiles()

	if len(mainModel.path_points) >= 2:
		curPathPointT = 0
		rotateToCurPathPoint()

	glutDisplayFunc(display)
	glutMouseFunc(onClick)
	glutMotionFunc(lambda x, y: arcBall.onMouseMove(x, y))
	glutReshapeFunc(onResize)
	glutKeyboardFunc(onKey)

	glutMainLoop()

if __name__ == "__main__":
	main()

