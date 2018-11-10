#########################################
# (c) 2011 Evan Mallory, TSM Games
#

# Input file format (text, which is output partly by
#    /blender/export_object.py and partly by
#    /blender/export_selected_vertices.py, or by
#    /python/makeModelFromSineMap.py):
# The first line is either "default", "spherical", or "pre_computed". If
#   it is "pre_computed", then there follow as many lines as there are
#   vertices, each containing two float values representing the texture
#   coordinates u and v of the ith vertex.
# <blank line>
# v1x v1y v1z n_v1x n_v1y n_y1z
#     .... v2 ....
#     .... vn ....
# <blank line>
# A line with the text "TRIANGLES" or "TRIANGLE_STRIP", specifying
#   what the following indices represent.
# <blank line>
# vertex_i1 vertex_i2
# vertex_i3 vertex_i4
# ... vertex_in ...
# <blank line>
# path_point_1_x path_point_1_y path_point_1_z path_point_1_type
# ... path_point_n_ ...
#
# Output file format (binary, native byte order):
# uint16: number of vertices
# For each vertex, 6 4-byte floats: x, y, z, nx, ny, nz
#           then 2 4-byte floats: tex_s, tex_t
# bool: true if indices represent TRIANGLE_STRIPS, false if TRIANGLES.
# uint32: number of vertex indices
# For each index, 1 uint16s: vertex_i
# uint16: number of path points
# For each path point, 3 4-byte floats: x, y, z
#                        then 1 uint16: type
# Finally, two uint8s indicating whether to clamp or repeat textures:
#    clamp_tex_s, clamp_tex_t

from utils import *

import math
import os
import re
from struct import pack
import sys

TEX_COORD_METHOD_DEFAULT = 'default'
TEX_COORD_METHOD_SPHERICAL = 'spherical'
TEX_COORD_METHOD_PRE_COMPUTED = 'pre_computed'
TEX_COORD_METHOD_IGNORE_Z = 'ignore_z'
TEX_COORD_METHOD_SINUSOIDAL = 'sinusoidal'
TEX_COORD_METHODS = [TEX_COORD_METHOD_DEFAULT,
					 TEX_COORD_METHOD_SPHERICAL,
					 TEX_COORD_METHOD_PRE_COMPUTED,
					 TEX_COORD_METHOD_IGNORE_Z,
					 TEX_COORD_METHOD_SINUSOIDAL]

VERTEX_INDICES_TYPES = ["TRIANGLE_STRIP", "TRIANGLES"]

class FileData():
	def __init__(self, vertices, tex_coord_method, tex_coords,
		vertex_indices_type, vertex_indices, path_points):
		self.vertices = vertices
		self.tex_coord_method = tex_coord_method
		self.tex_coords = tex_coords
		self.vertex_indices_type = vertex_indices_type
		self.vertex_indices = vertex_indices
		self.path_points = path_points


def readInput(file_in):
	MODE_TEX_COORD_METHOD = 0
	MODE_VERTICES = 1
	MODE_VERTEX_INDICES_TYPE = 2
	MODE_VERTEX_INDICES = 3
	MODE_PATH_POINTS = 4
	MODE_COUNT = 5 # Keep this up to date if a new mode is added

	mode = 0

	vertices = []
	vertex_indices_type = ""
	vertex_indices = []
	path_points = []
	tex_coord_method = None
	tex_coords = []

	for line in file_in:
		line = line.strip()
		if line == "":
			if mode == MODE_COUNT - 1:
				break;
			mode += 1
		else:
			if mode == MODE_TEX_COORD_METHOD and not tex_coord_method:
				tex_coord_method = line
				assert tex_coord_method in TEX_COORD_METHODS, \
					"Expected tex coord method in %s, found '%s'" % (
						TEX_COORD_METHODS, tex_coord_method)
			elif mode == MODE_VERTEX_INDICES_TYPE and not vertex_indices_type:
				vertex_indices_type = line
				assert vertex_indices_type in VERTEX_INDICES_TYPES, \
					"Expected vertex indices type in %s, found '%s'" % (
						VERTEX_INDICES_TYPES, vertex_indices_type)
			else:
				splits = re.split('[ \t]+', line)
				if mode == MODE_TEX_COORD_METHOD:
					assert tex_coord_method == TEX_COORD_METHOD_PRE_COMPUTED
					tex_coords.append([float(x) for x in splits[:2]])
				if mode == MODE_VERTICES:
					vertices.append([float(x) for x in splits[:6]])
				elif mode == MODE_VERTEX_INDICES:
					vertex_indices.extend([int(x) for x in splits])
				elif mode == MODE_PATH_POINTS:
					path_points.append([float(x) for x in splits[:3]] +
						[int(splits[3])])

	if tex_coord_method == TEX_COORD_METHOD_PRE_COMPUTED:
		assert len(tex_coords) == len(vertices), (
			"Read %d tex coords and %d vertices" % (
				len(tex_coords), len(vertices)))

	return FileData(vertices, tex_coord_method, tex_coords, vertex_indices_type, vertex_indices, path_points)

def computeTexCoords(data):
	tex_coord_method = data.tex_coord_method
	for i in range(0, len(data.vertices)):
		vert = data.vertices[i]

		if tex_coord_method == TEX_COORD_METHOD_DEFAULT:
			vert.append(0) # "u" coord

			radius = math.sqrt(math.pow(vert[0], 2) +
				math.pow(vert[1], 2) +
				math.pow(vert[2], 2))

			v = math.acos(vert[2] / radius) / math.pi
			vert.append((radius - .9) * 10)
		elif tex_coord_method == TEX_COORD_METHOD_SPHERICAL:
			(r, theta, phi) = xyzToSphericalCoords(vert[0:3])
			vert.append(theta / math.pi) # u
			vert.append((phi + math.pi) / (2 * math.pi)) # v
		elif tex_coord_method == TEX_COORD_METHOD_PRE_COMPUTED:
			vert.extend(data.tex_coords[i])
		elif tex_coord_method == TEX_COORD_METHOD_IGNORE_Z:
			vert.append(vert[0])
			vert.append(vert[1])
		elif tex_coord_method == TEX_COORD_METHOD_SINUSOIDAL:
			(x, y) = xyzToSinusoidal2d(vert[0:3], 2000, 1000)
			vert.append(x / 2000.0)
			# Textures have the origin in the bottom left, not top left.
			vert.append((1000 - y) / 1000.0)

def writeOutput(file_out, data):
	file_out.write(pack('@H', len(data.vertices)))
	for v in data.vertices:
		file_out.write(pack('@ffffffff', *v))

	file_out.write(pack('@B', data.vertex_indices_type == 'TRIANGLE_STRIP'))

	file_out.write(pack('@I', len(data.vertex_indices)))
	for f in data.vertex_indices:
		file_out.write(pack('@H', f))

	file_out.write(pack('@H', len(data.path_points)))
	for p in data.path_points:
		file_out.write(pack('@fffH', *p))

	if data.tex_coord_method == TEX_COORD_METHOD_PRE_COMPUTED:
		# By ESP, we know that we should clamp s and repeat t. We know this
		# because the same we wrote makeModelFromSineMap.py
		file_out.write(pack('@BB', 1, 0))
	else:
		file_out.write(pack('@BB', 1, 1))


def main():
	assert len(sys.argv) == 3, "Usage: %s <fin> <fout>" % sys.argv[0]

	file_in = open(sys.argv[1], 'rb')

	data = readInput(file_in)
	computeTexCoords(data)

	file_out = open(sys.argv[2], 'wb')
	writeOutput(file_out, data)

if __name__ == "__main__":
    main()
