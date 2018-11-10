"""
ArcBall.py -- Math utilities, vector, matrix types and
ArcBall quaternion rotation class
From http://nehe.gamedev.net/data/lessons/lesson.asp?lesson=48
"""

import numpy
import copy
from math import sqrt
from OpenGL.GLUT import GLUT_RIGHT_BUTTON, GLUT_LEFT_BUTTON, \
	GLUT_UP, GLUT_DOWN, glutPostRedisplay

# assuming IEEE-754(GLfloat), which i believe has max precision of 7 bits
Epsilon = 1.0e-5

class ArcBall:
	def __init__(self, NewWidth, NewHeight):
		self.m_StartVec = Vector3f()
		self.m_EndVec = Vector3f()
		self.m_AdjustWidth = 1.0
		self.m_AdjustHeight = 1.0
		self.setBounds(NewWidth, NewHeight)

		self.m_InitialRotation = Matrix3f()
		self.m_DragRotation = Matrix3f()
		self.m_Transform = Matrix4f()
		self.m_IsDragging = False

	def setBounds(self, NewWidth, NewHeight):
 		# Set new bounds
		assert(NewWidth > 1.0 and NewHeight > 1.0), "Invalid width or height for bounds."
		# Set adjustment factor for width/height
		self.m_AdjustWidth = 1.0 /((NewWidth - 1.0) * 0.5)
		self.m_AdjustHeight = 1.0 /((NewHeight - 1.0) * 0.5)

	def mapToSphere(self, NewPt):
		# Given a new window coordinate, will modify NewVec in place
		X = 0
		Y = 1
		Z = 2

		NewVec = Vector3f()
		# Copy paramter into temp point
		TempPt = copy.copy(NewPt)
		# Adjust point coords and scale down to range of [-1 ... 1]
		TempPt [X] =(NewPt [X] * self.m_AdjustWidth) - 1.0
		TempPt [Y] = 1.0 -(NewPt [Y] * self.m_AdjustHeight)
		# Compute the square of the length of the vector to the point from the center
		length = numpy.dot(TempPt, TempPt)
		# If the point is mapped outside of the sphere...(length > radius squared)
		if(length > 1.0):
			# Compute a normalizing factor(radius / sqrt(length))
			norm    = 1.0 / sqrt(length)

			# Return the "normalized" vector, a point on the sphere
			NewVec [X] = TempPt [X] * norm
			NewVec [Y] = TempPt [Y] * norm
			NewVec [Z] = 0.0
		else:			# Else it's on the inside
        	# Return a vector to a point mapped inside the sphere sqrt(radius squared - length)
			NewVec [X] = TempPt [X]
			NewVec [Y] = TempPt [Y]
			NewVec [Z] = sqrt(1.0 - length)

		return NewVec

	def getTransform(self):
		return self.m_Transform

	def setRotation(self, rotMatrix3):
		self.m_InitialRotation = numpy.array(rotMatrix3)
		self.m_Transform = Matrix4fSetRotationScaleFromMatrix3f(
			Matrix4f(), self.m_InitialRotation)

	def click(self, NewPt):
		# Mouse down(Point2fT
		self.m_StartVec = self.mapToSphere(NewPt)
		self.m_DragRotation = None

	def endDrag(self):
		if self.m_DragRotation != None:
			self.m_InitialRotation = numpy.array(self.m_DragRotation)

	def drag(self, NewPt):
		# Mouse drag, calculate rotation(Point2fQuat4fT)
		""" drag(Point2fT mouse_coord) -> new_quaternion_rotation_vec
		"""
		X = 0
		Y = 1
		Z = 2
		W = 3

		self.m_EndVec = self.mapToSphere(NewPt)

		# Compute the vector perpendicular to the begin and end vectors
		# Perp = Vector3f()
		Perp = numpy.cross(self.m_StartVec, self.m_EndVec)

		NewRot = Quat4f()
		# Compute the length of the perpendicular vector
		if(Vector3fLength(Perp) > Epsilon):		#    if its non-zero
			# We're ok, so return the perpendicular vector as the transform after all
			NewRot[X] = Perp[X]
			NewRot[Y] = Perp[Y]
			NewRot[Z] = Perp[Z]
			# In the quaternion values, w is cosine(theta / 2), where theta is rotation angle
			NewRot[W] = numpy.dot(self.m_StartVec, self.m_EndVec)
		else:		#                            if its zero
			# The begin and end vectors coincide, so return a quaternion of zero matrix(no rotation)
			NewRot[X] = NewRot[Y] = NewRot[Z] = NewRot[W] = 0.0

		self.m_DragRotation = Matrix3fSetRotationFromQuat4f(NewRot)
		self.m_DragRotation = numpy.dot(self.m_InitialRotation, self.m_DragRotation)
		drag_det = numpy.linalg.det(self.m_DragRotation)
		self.m_Transform = Matrix4fSetRotationScaleFromMatrix3f(
			Matrix4f(), self.m_DragRotation)


	# ##################### GL event utilities ##############################

	def onClick(self, button, button_state, x, y):
		self.m_IsDragging = False
		if (button == GLUT_RIGHT_BUTTON and button_state == GLUT_UP):
			# Right button click
			self.setRotation(Matrix3f())
			glutPostRedisplay()
		elif (button == GLUT_LEFT_BUTTON and button_state == GLUT_UP):
			# Left button clicked down
			self.m_IsDragging = False
			self.endDrag()
		elif (button == GLUT_LEFT_BUTTON and button_state == GLUT_DOWN):
			# Left button clicked down
			self.m_IsDragging = True
			self.click(Point2f(x, y))

	def onMouseMove(self, x, y):
		if self.m_IsDragging:
			mouse_pt = Point2f(x, y)
			self.drag(mouse_pt)
			glutPostRedisplay()
		return

# ##################### Math utilities#######################################


def Matrix4f():
	return numpy.identity(4, 'f')

def Matrix3f():
	return numpy.identity(3, 'f')

def Quat4f():
	return numpy.zeros(4, 'f')

def Vector3f():
	return numpy.zeros(3, 'f')

def Point2f(x = 0.0, y = 0.0):
	pt = numpy.zeros(2, 'f')
	pt [0] = x
	pt [1] = y
	return pt

def Vector3fLength(u):
	return sqrt(numpy.dot(u, u))

def Matrix4fSetRotationScaleFromMatrix3f(NewObj, three_by_three_matrix):
	# Modifies NewObj in-place by replacing its upper 3x3 portion from the
	# passed in 3x3 matrix.
	# NewObj = Matrix4f()
	NewObj [0:3,0:3] = three_by_three_matrix
	return NewObj

def Matrix3fSetRotationFromQuat4f(q1):
	# Converts the H quaternion q1 into a new equivalent 3x3 rotation matrix.
	X = 0
	Y = 1
	Z = 2
	W = 3

	NewObj = Matrix3f()
	n = numpy.dot(q1, q1)
	s = 0.0
	if(n > 0.0):
		s = 2.0 / n
	xs = q1 [X] * s;  ys = q1 [Y] * s;  zs = q1 [Z] * s;
	wx = q1 [W] * xs; wy = q1 [W] * ys; wz = q1 [W] * zs;
	xx = q1 [X] * xs; xy = q1 [X] * ys; xz = q1 [X] * zs;
	yy = q1 [Y] * ys; yz = q1 [Y] * zs; zz = q1 [Z] * zs;
	# This math all comes about by way of algebra, complex math, and trig identities.
	# See Lengyel pages 88-92
	NewObj [X][X] = 1.0 -(yy + zz);	NewObj [Y][X] = xy - wz; 			NewObj [Z][X] = xz + wy;
	NewObj [X][Y] =       xy + wz;	NewObj [Y][Y] = 1.0 -(xx + zz);		NewObj [Z][Y] = yz - wx;
	NewObj [X][Z] =       xz - wy;	NewObj [Y][Z] = yz + wx;          	NewObj [Z][Z] = 1.0 -(xx + yy);

	return NewObj

