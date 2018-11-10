#########################################
# (c) 2011 Evan Mallory, TSM Games
#

import sys
from utils import *
from math import sqrt

# This creates a square that has an edge that is slightly lower than 0 z
# and everything but the edge is slightly higher than 0 z. It needs
# to fit over a piece of the (spherical) path without showing any gaps
# and without being covered by the path. The part of the square
# that is above 0 z is in the shape of a circle.

def main ():
	assert len(sys.argv) == 2,\
		"Usage: %s <output_file>" % sys.argv[0]

	num_divisions = 8

	fout = open(sys.argv[1], 'w')
	fout.write("pre_computed\n")
	for y in drange(0, 1.001, 1.0 / num_divisions):
		for x in drange(0, 1.001, 1.0 / num_divisions):
			fout.write("%f %f\n" % (y, x))
	fout.write("\n")

	triangles = []
	vertex_i = 0

	scale = .14

	for y in drange(-.5, .501, 1.0 / num_divisions):
		for x in drange(-.5, .501, 1.0 / num_divisions):
			z = .005
			if sqrt(x*x + y*y) > .5:
				z = -.01
			fout.write("%f %f %f %f %f %f\n" % (x*scale, y*scale, z, 0, 0, 1))
			if y > -.5:
				if x > -.5:
					triangles.extend([vertex_i, vertex_i - 1,
									  vertex_i - (num_divisions + 1)])
				if x < .5:
					triangles.extend([vertex_i,
									  vertex_i - (num_divisions + 1),
									  vertex_i - num_divisions])
			vertex_i += 1

	assert len(triangles)/3 == num_divisions*num_divisions * 2

	fout.write("\n")
	for i in range(0, len(triangles) / 3):
		fout.write("%d %d %d\n" % (triangles[3*i], triangles[3*i+1], triangles[3*i+2]))


if __name__ == "__main__":
	main()
