#!/usr/bin/env python2

#########################################
# (c) 2011 Evan Mallory, TSM Games
# (c) 2018 James Dowdell, TSM Games
#########################################


from struct import unpack
from OpenGL.GL import *

VERTEX_POINT_SIZE = 32
PATH_POINT_SIZE = 3 * 4 + 2  # 3 4-byte floats, 1 uint16

def drawModel(m, texture=None):
    glEnableClientState(GL_VERTEX_ARRAY)
    glEnableClientState(GL_NORMAL_ARRAY)
    if texture:
        glEnable(GL_TEXTURE_2D)
        glEnableClientState(GL_TEXTURE_COORD_ARRAY)
        glBindTexture(GL_TEXTURE_2D, texture)

    glVertexPointer(3, GL_FLOAT, VERTEX_POINT_SIZE, m.vertices)
    glNormalPointer(GL_FLOAT, VERTEX_POINT_SIZE, m.vertices[12:])
    if texture:
        glTexCoordPointer(2, GL_FLOAT, VERTEX_POINT_SIZE, m.vertices[24:])

    if m.use_triangle_strip:
        glDrawElements(GL_TRIANGLE_STRIP, len(m.vertex_indices) / 2, GL_UNSIGNED_SHORT, m.vertex_indices)
    else:
        glDrawElements(GL_TRIANGLES, len(m.vertex_indices) / 2, GL_UNSIGNED_SHORT, m.vertex_indices)

    glDisableClientState(GL_VERTEX_ARRAY)
    glDisableClientState(GL_NORMAL_ARRAY)

    glDisableClientState(GL_TEXTURE_COORD_ARRAY)
    glDisable(GL_TEXTURE_2D)


class ModelData():
    def __init__(self, vertices, vertex_indices, use_triangle_strip, path_points, clampS, clampT):
        # vertices and vertex_indices are binary data
        self.vertices = vertices
        self.vertex_indices = vertex_indices
        self.use_triangle_strip = use_triangle_strip
        # path_points is an array [[float, float, float, short], ... ]
        self.path_points = path_points
        self.clampS = clampS
        self.clampT = clampT

def loadModel(filename):
    file_in = open(filename, 'rb')
    num_vertices = unpack('@H', file_in.read(2))[0]
    vertices_normals_tex = file_in.read(num_vertices * VERTEX_POINT_SIZE)

    for i in range(0, min(10, num_vertices)):
        unpacked = unpack('@ffffffff',
            vertices_normals_tex[i*VERTEX_POINT_SIZE:(i+1)*VERTEX_POINT_SIZE])
        print ("%6.3f "*6 +"  " + "%.2f "*2) % unpacked

    use_triangle_strip = unpack('@b', file_in.read(1))[0]
    num_vertex_indices = unpack('@I', file_in.read(4))[0]
    vertex_indices = file_in.read(num_vertex_indices * 2)
    for i in range(0, min(10, num_vertex_indices)):
        unpacked = unpack('@H', vertex_indices[i*2:(i+1)*2])
        print unpacked[0]

    num_path_points = unpack('@H', file_in.read(2))[0]
    num_waypoints = 0
    path_point_data = file_in.read(num_path_points * PATH_POINT_SIZE)
    path_points = []
    for i in range(0, num_path_points):
        unpacked = unpack('@fffH',
            path_point_data[i*PATH_POINT_SIZE:(i+1)*PATH_POINT_SIZE])
        path_points.append(unpacked)
        if i < 10:
            print ("%5.2f " * 3 + "%2d") % unpacked
        if unpacked[3] == 1:
            num_waypoints += 1

    clampTextureData = file_in.read(2)
    (clampS, clampT) = unpack('@BB', clampTextureData)

    print "Read %d vertices, %d vertex indices (%s), and %d path points (%d waypoints)" % (
        num_vertices, num_vertex_indices,
        use_triangle_strip and "TRIANGLE_STRIP" or "TRIANGLES",
        num_path_points, num_waypoints)

    return ModelData(vertices_normals_tex, vertex_indices, use_triangle_strip, path_points, clampS, clampT)


###############################################################
#                     OpenGL debugging
###############################################################

def getLightParams():
    parameters = [
        (GL_AMBIENT, "GL_AMBIENT"),
        (GL_DIFFUSE, "GL_DIFFUSE"),
        (GL_SPECULAR, "GL_SPECULAR"),
        (GL_POSITION, "GL_POSITION"),
        (GL_SPOT_DIRECTION, "GL_SPOT_DIRECTION"),
        (GL_SPOT_EXPONENT,"GL_SPOT_EXPONENT"),
        (GL_SPOT_CUTOFF,"GL_SPOT_CUTOFF"),
        (GL_CONSTANT_ATTENUATION,"GL_CONSTANT_ATTENUATION"),
        (GL_LINEAR_ATTENUATION , "GL_LINEAR_ATTENUATION"),
        (GL_QUADRATIC_ATTENUATION, "GL_QUADRATIC_ATTENUATION"),
    ]

    for param, name in parameters:
        print name, glGetLightfv(GL_LIGHT0, param ) # now requires fully-specified name...

def getMaterialParams():
    parameters = [
        (GL_AMBIENT, "GL_AMBIENT"),
        (GL_DIFFUSE, "GL_DIFFUSE"),
        (GL_SPECULAR, "GL_SPECULAR"),
        (GL_EMISSION, "GL_EMISSION"),
        (GL_SHININESS, "GL_SHININESS"),
        (GL_COLOR_INDEXES, "GL_COLOR_INDEXES"),
    ]

    print
    for param, name in parameters:
        print name, glGetMaterialfv(GL_FRONT, param)


# Excluded:
# GL_CLIP_PLANE i
# GL_LIGHT i
def getGetParams():
    parameters = [
        (GL_ACCUM_ALPHA_BITS, "GL_ACCUM_ALPHA_BITS"),
        (GL_ACCUM_BLUE_BITS, "GL_ACCUM_BLUE_BITS"),
        (GL_ACCUM_CLEAR_VALUE, "GL_ACCUM_CLEAR_VALUE"),
        (GL_ACCUM_GREEN_BITS, "GL_ACCUM_GREEN_BITS"),
        (GL_ACCUM_RED_BITS, "GL_ACCUM_RED_BITS"),
        (GL_ACTIVE_TEXTURE_ARB, "GL_ACTIVE_TEXTURE_ARB"),
        (GL_ALIASED_POINT_SIZE_RANGE, "GL_ALIASED_POINT_SIZE_RANGE"),
        (GL_ALIASED_LINE_WIDTH_RANGE, "GL_ALIASED_LINE_WIDTH_RANGE"),
        (GL_ALPHA_BIAS, "GL_ALPHA_BIAS"),
        (GL_ALPHA_BITS, "GL_ALPHA_BITS"),
        (GL_ALPHA_SCALE, "GL_ALPHA_SCALE"),
        (GL_ALPHA_TEST, "GL_ALPHA_TEST"),
        (GL_ALPHA_TEST_FUNC, "GL_ALPHA_TEST_FUNC"),
        (GL_ALPHA_TEST_REF, "GL_ALPHA_TEST_REF"),
        (GL_ATTRIB_STACK_DEPTH, "GL_ATTRIB_STACK_DEPTH"),
        (GL_AUTO_NORMAL, "GL_AUTO_NORMAL"),
        (GL_AUX_BUFFERS, "GL_AUX_BUFFERS"),
        (GL_BLEND, "GL_BLEND"),
        (GL_BLEND_COLOR, "GL_BLEND_COLOR"),
        (GL_BLEND_DST, "GL_BLEND_DST"),
        (GL_BLEND_EQUATION, "GL_BLEND_EQUATION"),
        (GL_BLEND_SRC, "GL_BLEND_SRC"),
        (GL_BLUE_BIAS, "GL_BLUE_BIAS"),
        (GL_BLUE_BITS, "GL_BLUE_BITS"),
        (GL_BLUE_SCALE, "GL_BLUE_SCALE"),
        (GL_CLIENT_ACTIVE_TEXTURE_ARB, "GL_CLIENT_ACTIVE_TEXTURE_ARB"),
        (GL_CLIENT_ATTRIB_STACK_DEPTH, "GL_CLIENT_ATTRIB_STACK_DEPTH"),
        (GL_COLOR_ARRAY, "GL_COLOR_ARRAY"),
        (GL_COLOR_ARRAY_SIZE, "GL_COLOR_ARRAY_SIZE"),
        (GL_COLOR_ARRAY_STRIDE, "GL_COLOR_ARRAY_STRIDE"),
        (GL_COLOR_ARRAY_TYPE, "GL_COLOR_ARRAY_TYPE"),
        (GL_COLOR_CLEAR_VALUE, "GL_COLOR_CLEAR_VALUE"),
        (GL_COLOR_LOGIC_OP, "GL_COLOR_LOGIC_OP"),
        (GL_COLOR_MATERIAL, "GL_COLOR_MATERIAL"),
        (GL_COLOR_MATERIAL_FACE, "GL_COLOR_MATERIAL_FACE"),
        (GL_COLOR_MATERIAL_PARAMETER, "GL_COLOR_MATERIAL_PARAMETER"),
        (GL_COLOR_MATRIX, "GL_COLOR_MATRIX"),
        (GL_COLOR_MATRIX_STACK_DEPTH, "GL_COLOR_MATRIX_STACK_DEPTH"),
        (GL_COLOR_TABLE, "GL_COLOR_TABLE"),
        (GL_COLOR_WRITEMASK, "GL_COLOR_WRITEMASK"),
        (GL_CONVOLUTION_1D, "GL_CONVOLUTION_1D"),
        (GL_CONVOLUTION_2D, "GL_CONVOLUTION_2D"),
        (GL_CULL_FACE, "GL_CULL_FACE"),
        (GL_CULL_FACE_MODE, "GL_CULL_FACE_MODE"),
        (GL_CURRENT_COLOR, "GL_CURRENT_COLOR"),
        (GL_CURRENT_INDEX, "GL_CURRENT_INDEX"),
        (GL_CURRENT_NORMAL, "GL_CURRENT_NORMAL"),
        (GL_CURRENT_RASTER_COLOR, "GL_CURRENT_RASTER_COLOR"),
        (GL_CURRENT_RASTER_DISTANCE, "GL_CURRENT_RASTER_DISTANCE"),
        (GL_CURRENT_RASTER_INDEX, "GL_CURRENT_RASTER_INDEX"),
        (GL_CURRENT_RASTER_POSITION, "GL_CURRENT_RASTER_POSITION"),
        (GL_CURRENT_RASTER_POSITION_VALID, "GL_CURRENT_RASTER_POSITION_VALID"),
        (GL_CURRENT_RASTER_TEXTURE_COORDS, "GL_CURRENT_RASTER_TEXTURE_COORDS"),
        (GL_CURRENT_TEXTURE_COORDS, "GL_CURRENT_TEXTURE_COORDS"),
        (GL_DEPTH_BIAS, "GL_DEPTH_BIAS"),
        (GL_DEPTH_BITS, "GL_DEPTH_BITS"),
        (GL_DEPTH_CLEAR_VALUE, "GL_DEPTH_CLEAR_VALUE"),
        (GL_DEPTH_FUNC, "GL_DEPTH_FUNC"),
        (GL_DEPTH_RANGE, "GL_DEPTH_RANGE"),
        (GL_DEPTH_SCALE, "GL_DEPTH_SCALE"),
        (GL_DEPTH_TEST, "GL_DEPTH_TEST"),
        (GL_DEPTH_WRITEMASK, "GL_DEPTH_WRITEMASK"),
        (GL_DITHER, "GL_DITHER"),
        (GL_DOUBLEBUFFER, "GL_DOUBLEBUFFER"),
        (GL_DRAW_BUFFER, "GL_DRAW_BUFFER"),
        (GL_EDGE_FLAG, "GL_EDGE_FLAG"),
        (GL_EDGE_FLAG_ARRAY, "GL_EDGE_FLAG_ARRAY"),
        (GL_EDGE_FLAG_ARRAY_STRIDE, "GL_EDGE_FLAG_ARRAY_STRIDE"),
        (GL_FEEDBACK_BUFFER_SIZE, "GL_FEEDBACK_BUFFER_SIZE"),
        (GL_FEEDBACK_BUFFER_TYPE, "GL_FEEDBACK_BUFFER_TYPE"),
        (GL_FOG, "GL_FOG"),
        (GL_FOG_COLOR, "GL_FOG_COLOR"),
        (GL_FOG_DENSITY, "GL_FOG_DENSITY"),
        (GL_FOG_END, "GL_FOG_END"),
        (GL_FOG_HINT, "GL_FOG_HINT"),
        (GL_FOG_INDEX, "GL_FOG_INDEX"),
        (GL_FOG_MODE, "GL_FOG_MODE"),
        (GL_FOG_START, "GL_FOG_START"),
        (GL_FRONT_FACE, "GL_FRONT_FACE"),
        (GL_GREEN_BIAS, "GL_GREEN_BIAS"),
        (GL_GREEN_BITS, "GL_GREEN_BITS"),
        (GL_GREEN_SCALE, "GL_GREEN_SCALE"),
        (GL_HISTOGRAM, "GL_HISTOGRAM"),
        (GL_INDEX_ARRAY, "GL_INDEX_ARRAY"),
        (GL_INDEX_ARRAY_STRIDE, "GL_INDEX_ARRAY_STRIDE"),
        (GL_INDEX_ARRAY_TYPE, "GL_INDEX_ARRAY_TYPE"),
        (GL_INDEX_BITS, "GL_INDEX_BITS"),
        (GL_INDEX_CLEAR_VALUE, "GL_INDEX_CLEAR_VALUE"),
        (GL_INDEX_LOGIC_OP, "GL_INDEX_LOGIC_OP"),
        (GL_INDEX_MODE, "GL_INDEX_MODE"),
        (GL_INDEX_OFFSET, "GL_INDEX_OFFSET"),
        (GL_INDEX_SHIFT, "GL_INDEX_SHIFT"),
        (GL_INDEX_WRITEMASK, "GL_INDEX_WRITEMASK"),
        (GL_LIGHTING, "GL_LIGHTING"),
        (GL_LIGHT_MODEL_AMBIENT, "GL_LIGHT_MODEL_AMBIENT"),
        (GL_LIGHT_MODEL_COLOR_CONTROL, "GL_LIGHT_MODEL_COLOR_CONTROL"),
        (GL_LIGHT_MODEL_LOCAL_VIEWER, "GL_LIGHT_MODEL_LOCAL_VIEWER"),
        (GL_LIGHT_MODEL_TWO_SIDE, "GL_LIGHT_MODEL_TWO_SIDE"),
        (GL_LINE_SMOOTH, "GL_LINE_SMOOTH"),
        (GL_LINE_SMOOTH_HINT, "GL_LINE_SMOOTH_HINT"),
        (GL_LINE_STIPPLE, "GL_LINE_STIPPLE"),
        (GL_LINE_STIPPLE_PATTERN, "GL_LINE_STIPPLE_PATTERN"),
        (GL_LINE_STIPPLE_REPEAT, "GL_LINE_STIPPLE_REPEAT"),
        (GL_LINE_WIDTH, "GL_LINE_WIDTH"),
        (GL_LINE_WIDTH_GRANULARITY, "GL_LINE_WIDTH_GRANULARITY"),
        (GL_LINE_WIDTH_RANGE, "GL_LINE_WIDTH_RANGE"),
        (GL_LIST_BASE, "GL_LIST_BASE"),
        (GL_LIST_INDEX, "GL_LIST_INDEX"),
        (GL_LIST_MODE, "GL_LIST_MODE"),
        (GL_LOGIC_OP_MODE, "GL_LOGIC_OP_MODE"),
        (GL_MAP1_COLOR_4, "GL_MAP1_COLOR_4"),
        (GL_MAP1_GRID_DOMAIN, "GL_MAP1_GRID_DOMAIN"),
        (GL_MAP1_GRID_SEGMENTS, "GL_MAP1_GRID_SEGMENTS"),
        (GL_MAP1_INDEX, "GL_MAP1_INDEX"),
        (GL_MAP1_NORMAL, "GL_MAP1_NORMAL"),
        (GL_MAP1_TEXTURE_COORD_1, "GL_MAP1_TEXTURE_COORD_1"),
        (GL_MAP1_TEXTURE_COORD_2, "GL_MAP1_TEXTURE_COORD_2"),
        (GL_MAP1_TEXTURE_COORD_3, "GL_MAP1_TEXTURE_COORD_3"),
        (GL_MAP1_TEXTURE_COORD_4, "GL_MAP1_TEXTURE_COORD_4"),
        (GL_MAP1_VERTEX_3, "GL_MAP1_VERTEX_3"),
        (GL_MAP1_VERTEX_4, "GL_MAP1_VERTEX_4"),
        (GL_MAP2_COLOR_4, "GL_MAP2_COLOR_4"),
        (GL_MAP2_GRID_DOMAIN, "GL_MAP2_GRID_DOMAIN"),
        (GL_MAP2_GRID_SEGMENTS, "GL_MAP2_GRID_SEGMENTS"),
        (GL_MAP2_INDEX, "GL_MAP2_INDEX"),
        (GL_MAP2_NORMAL, "GL_MAP2_NORMAL"),
        (GL_MAP2_TEXTURE_COORD_1, "GL_MAP2_TEXTURE_COORD_1"),
        (GL_MAP2_TEXTURE_COORD_2, "GL_MAP2_TEXTURE_COORD_2"),
        (GL_MAP2_TEXTURE_COORD_3, "GL_MAP2_TEXTURE_COORD_3"),
        (GL_MAP2_TEXTURE_COORD_4, "GL_MAP2_TEXTURE_COORD_4"),
        (GL_MAP2_VERTEX_3, "GL_MAP2_VERTEX_3"),
        (GL_MAP2_VERTEX_4, "GL_MAP2_VERTEX_4"),
        (GL_MAP_COLOR, "GL_MAP_COLOR"),
        (GL_MAP_STENCIL, "GL_MAP_STENCIL"),
        (GL_MATRIX_MODE, "GL_MATRIX_MODE"),
        (GL_MAX_3D_TEXTURE_SIZE, "GL_MAX_3D_TEXTURE_SIZE"),
        (GL_MAX_CLIENT_ATTRIB_STACK_DEPTH, "GL_MAX_CLIENT_ATTRIB_STACK_DEPTH"),
        (GL_MAX_ATTRIB_STACK_DEPTH, "GL_MAX_ATTRIB_STACK_DEPTH"),
        (GL_MAX_CLIP_PLANES, "GL_MAX_CLIP_PLANES"),
        (GL_MAX_COLOR_MATRIX_STACK_DEPTH, "GL_MAX_COLOR_MATRIX_STACK_DEPTH"),
        (GL_MAX_ELEMENTS_INDICES, "GL_MAX_ELEMENTS_INDICES"),
        (GL_MAX_ELEMENTS_VERTICES, "GL_MAX_ELEMENTS_VERTICES"),
        (GL_MAX_EVAL_ORDER, "GL_MAX_EVAL_ORDER"),
        (GL_MAX_LIGHTS, "GL_MAX_LIGHTS"),
        (GL_MAX_LIST_NESTING, "GL_MAX_LIST_NESTING"),
        (GL_MAX_MODELVIEW_STACK_DEPTH, "GL_MAX_MODELVIEW_STACK_DEPTH"),
        (GL_MAX_NAME_STACK_DEPTH, "GL_MAX_NAME_STACK_DEPTH"),
        (GL_MAX_PIXEL_MAP_TABLE, "GL_MAX_PIXEL_MAP_TABLE"),
        (GL_MAX_PROJECTION_STACK_DEPTH, "GL_MAX_PROJECTION_STACK_DEPTH"),
        (GL_MAX_TEXTURE_SIZE, "GL_MAX_TEXTURE_SIZE"),
        (GL_MAX_TEXTURE_STACK_DEPTH, "GL_MAX_TEXTURE_STACK_DEPTH"),
        (GL_MAX_TEXTURE_UNITS_ARB, "GL_MAX_TEXTURE_UNITS_ARB"),
        (GL_MAX_VIEWPORT_DIMS, "GL_MAX_VIEWPORT_DIMS"),
        (GL_MINMAX, "GL_MINMAX"),
        (GL_MODELVIEW_MATRIX, "GL_MODELVIEW_MATRIX"),
        (GL_MODELVIEW_STACK_DEPTH, "GL_MODELVIEW_STACK_DEPTH"),
        (GL_NAME_STACK_DEPTH, "GL_NAME_STACK_DEPTH"),
        (GL_NORMAL_ARRAY, "GL_NORMAL_ARRAY"),
        (GL_NORMAL_ARRAY_STRIDE, "GL_NORMAL_ARRAY_STRIDE"),
        (GL_NORMAL_ARRAY_TYPE, "GL_NORMAL_ARRAY_TYPE"),
        (GL_NORMALIZE, "GL_NORMALIZE"),
        (GL_PACK_ALIGNMENT, "GL_PACK_ALIGNMENT"),
        (GL_PACK_IMAGE_HEIGHT, "GL_PACK_IMAGE_HEIGHT"),
        (GL_PACK_LSB_FIRST, "GL_PACK_LSB_FIRST"),
        (GL_PACK_ROW_LENGTH, "GL_PACK_ROW_LENGTH"),
        (GL_PACK_SKIP_IMAGES, "GL_PACK_SKIP_IMAGES"),
        (GL_PACK_SKIP_PIXELS, "GL_PACK_SKIP_PIXELS"),
        (GL_PACK_SKIP_ROWS, "GL_PACK_SKIP_ROWS"),
        (GL_PACK_SWAP_BYTES, "GL_PACK_SWAP_BYTES"),
        (GL_PERSPECTIVE_CORRECTION_HINT, "GL_PERSPECTIVE_CORRECTION_HINT"),
        (GL_PIXEL_MAP_A_TO_A_SIZE, "GL_PIXEL_MAP_A_TO_A_SIZE"),
        (GL_PIXEL_MAP_B_TO_B_SIZE, "GL_PIXEL_MAP_B_TO_B_SIZE"),
        (GL_PIXEL_MAP_G_TO_G_SIZE, "GL_PIXEL_MAP_G_TO_G_SIZE"),
        (GL_PIXEL_MAP_I_TO_A_SIZE, "GL_PIXEL_MAP_I_TO_A_SIZE"),
        (GL_PIXEL_MAP_I_TO_B_SIZE, "GL_PIXEL_MAP_I_TO_B_SIZE"),
        (GL_PIXEL_MAP_I_TO_G_SIZE, "GL_PIXEL_MAP_I_TO_G_SIZE"),
        (GL_PIXEL_MAP_I_TO_I_SIZE, "GL_PIXEL_MAP_I_TO_I_SIZE"),
        (GL_PIXEL_MAP_I_TO_R_SIZE, "GL_PIXEL_MAP_I_TO_R_SIZE"),
        (GL_PIXEL_MAP_R_TO_R_SIZE, "GL_PIXEL_MAP_R_TO_R_SIZE"),
        (GL_PIXEL_MAP_S_TO_S_SIZE, "GL_PIXEL_MAP_S_TO_S_SIZE"),
        (GL_POINT_SIZE, "GL_POINT_SIZE"),
        (GL_POINT_SIZE_GRANULARITY, "GL_POINT_SIZE_GRANULARITY"),
        (GL_POINT_SIZE_RANGE, "GL_POINT_SIZE_RANGE"),
        (GL_POINT_SMOOTH, "GL_POINT_SMOOTH"),
        (GL_POINT_SMOOTH_HINT, "GL_POINT_SMOOTH_HINT"),
        (GL_POLYGON_MODE, "GL_POLYGON_MODE"),
        (GL_POLYGON_OFFSET_FACTOR, "GL_POLYGON_OFFSET_FACTOR"),
        (GL_POLYGON_OFFSET_UNITS, "GL_POLYGON_OFFSET_UNITS"),
        (GL_POLYGON_OFFSET_FILL, "GL_POLYGON_OFFSET_FILL"),
        (GL_POLYGON_OFFSET_LINE, "GL_POLYGON_OFFSET_LINE"),
        (GL_POLYGON_OFFSET_POINT, "GL_POLYGON_OFFSET_POINT"),
        (GL_POLYGON_SMOOTH, "GL_POLYGON_SMOOTH"),
        (GL_POLYGON_SMOOTH_HINT, "GL_POLYGON_SMOOTH_HINT"),
        (GL_POLYGON_STIPPLE, "GL_POLYGON_STIPPLE"),
        (GL_POST_COLOR_MATRIX_COLOR_TABLE, "GL_POST_COLOR_MATRIX_COLOR_TABLE"),
        (GL_POST_COLOR_MATRIX_RED_BIAS, "GL_POST_COLOR_MATRIX_RED_BIAS"),
        (GL_POST_COLOR_MATRIX_GREEN_BIAS, "GL_POST_COLOR_MATRIX_GREEN_BIAS"),
        (GL_POST_COLOR_MATRIX_BLUE_BIAS, "GL_POST_COLOR_MATRIX_BLUE_BIAS"),
        (GL_POST_COLOR_MATRIX_ALPHA_BIAS, "GL_POST_COLOR_MATRIX_ALPHA_BIAS"),
        (GL_POST_COLOR_MATRIX_RED_SCALE, "GL_POST_COLOR_MATRIX_RED_SCALE"),
        (GL_POST_COLOR_MATRIX_GREEN_SCALE, "GL_POST_COLOR_MATRIX_GREEN_SCALE"),
        (GL_POST_COLOR_MATRIX_BLUE_SCALE, "GL_POST_COLOR_MATRIX_BLUE_SCALE"),
        (GL_POST_COLOR_MATRIX_ALPHA_SCALE, "GL_POST_COLOR_MATRIX_ALPHA_SCALE"),
        (GL_POST_CONVOLUTION_COLOR_TABLE, "GL_POST_CONVOLUTION_COLOR_TABLE"),
        (GL_POST_CONVOLUTION_RED_BIAS, "GL_POST_CONVOLUTION_RED_BIAS"),
        (GL_POST_CONVOLUTION_GREEN_BIAS, "GL_POST_CONVOLUTION_GREEN_BIAS"),
        (GL_POST_CONVOLUTION_BLUE_BIAS, "GL_POST_CONVOLUTION_BLUE_BIAS"),
        (GL_POST_CONVOLUTION_ALPHA_BIAS, "GL_POST_CONVOLUTION_ALPHA_BIAS"),
        (GL_POST_CONVOLUTION_RED_SCALE, "GL_POST_CONVOLUTION_RED_SCALE"),
        (GL_POST_CONVOLUTION_GREEN_SCALE, "GL_POST_CONVOLUTION_GREEN_SCALE"),
        (GL_POST_CONVOLUTION_BLUE_SCALE, "GL_POST_CONVOLUTION_BLUE_SCALE"),
        (GL_POST_CONVOLUTION_ALPHA_SCALE, "GL_POST_CONVOLUTION_ALPHA_SCALE"),
        (GL_PROJECTION_MATRIX, "GL_PROJECTION_MATRIX"),
        (GL_PROJECTION_STACK_DEPTH, "GL_PROJECTION_STACK_DEPTH"),
        (GL_READ_BUFFER, "GL_READ_BUFFER"),
        (GL_RED_BIAS, "GL_RED_BIAS"),
        (GL_RED_BITS, "GL_RED_BITS"),
        (GL_RED_SCALE, "GL_RED_SCALE"),
        (GL_RENDER_MODE, "GL_RENDER_MODE"),
        (GL_RESCALE_NORMAL, "GL_RESCALE_NORMAL"),
        (GL_RGBA_MODE, "GL_RGBA_MODE"),
        (GL_SCISSOR_BOX, "GL_SCISSOR_BOX"),
        (GL_SCISSOR_TEST, "GL_SCISSOR_TEST"),
        (GL_SELECTION_BUFFER_SIZE, "GL_SELECTION_BUFFER_SIZE"),
        (GL_SEPARABLE_2D, "GL_SEPARABLE_2D"),
        (GL_SHADE_MODEL, "GL_SHADE_MODEL"),
        (GL_SMOOTH_LINE_WIDTH_RANGE, "GL_SMOOTH_LINE_WIDTH_RANGE"),
        (GL_SMOOTH_LINE_WIDTH_GRANULARITY, "GL_SMOOTH_LINE_WIDTH_GRANULARITY"),
        (GL_SMOOTH_POINT_SIZE_RANGE, "GL_SMOOTH_POINT_SIZE_RANGE"),
        (GL_SMOOTH_POINT_SIZE_GRANULARITY, "GL_SMOOTH_POINT_SIZE_GRANULARITY"),
        (GL_STENCIL_BITS, "GL_STENCIL_BITS"),
        (GL_STENCIL_CLEAR_VALUE, "GL_STENCIL_CLEAR_VALUE"),
        (GL_STENCIL_FAIL, "GL_STENCIL_FAIL"),
        (GL_STENCIL_FUNC, "GL_STENCIL_FUNC"),
        (GL_STENCIL_PASS_DEPTH_FAIL, "GL_STENCIL_PASS_DEPTH_FAIL"),
        (GL_STENCIL_PASS_DEPTH_PASS, "GL_STENCIL_PASS_DEPTH_PASS"),
        (GL_STENCIL_REF, "GL_STENCIL_REF"),
        (GL_STENCIL_TEST, "GL_STENCIL_TEST"),
        (GL_STENCIL_VALUE_MASK, "GL_STENCIL_VALUE_MASK"),
        (GL_STENCIL_WRITEMASK, "GL_STENCIL_WRITEMASK"),
        (GL_STEREO, "GL_STEREO"),
        (GL_SUBPIXEL_BITS, "GL_SUBPIXEL_BITS"),
        (GL_TEXTURE_1D, "GL_TEXTURE_1D"),
        (GL_TEXTURE_BINDING_1D, "GL_TEXTURE_BINDING_1D"),
        (GL_TEXTURE_2D, "GL_TEXTURE_2D"),
        (GL_TEXTURE_BINDING_2D, "GL_TEXTURE_BINDING_2D"),
        (GL_TEXTURE_3D, "GL_TEXTURE_3D"),
        (GL_TEXTURE_BINDING_3D, "GL_TEXTURE_BINDING_3D"),
        (GL_TEXTURE_COORD_ARRAY, "GL_TEXTURE_COORD_ARRAY"),
        (GL_TEXTURE_COORD_ARRAY_SIZE, "GL_TEXTURE_COORD_ARRAY_SIZE"),
        (GL_TEXTURE_COORD_ARRAY_STRIDE, "GL_TEXTURE_COORD_ARRAY_STRIDE"),
        (GL_TEXTURE_COORD_ARRAY_TYPE, "GL_TEXTURE_COORD_ARRAY_TYPE"),
        (GL_TEXTURE_GEN_Q, "GL_TEXTURE_GEN_Q"),
        (GL_TEXTURE_GEN_R, "GL_TEXTURE_GEN_R"),
        (GL_TEXTURE_GEN_S, "GL_TEXTURE_GEN_S"),
        (GL_TEXTURE_GEN_T, "GL_TEXTURE_GEN_T"),
        (GL_TEXTURE_MATRIX, "GL_TEXTURE_MATRIX"),
        (GL_TEXTURE_STACK_DEPTH, "GL_TEXTURE_STACK_DEPTH"),
        (GL_UNPACK_ALIGNMENT, "GL_UNPACK_ALIGNMENT"),
        (GL_UNPACK_IMAGE_HEIGHT, "GL_UNPACK_IMAGE_HEIGHT"),
        (GL_UNPACK_LSB_FIRST, "GL_UNPACK_LSB_FIRST"),
        (GL_UNPACK_ROW_LENGTH, "GL_UNPACK_ROW_LENGTH"),
        (GL_UNPACK_SKIP_IMAGES, "GL_UNPACK_SKIP_IMAGES"),
        (GL_UNPACK_SKIP_PIXELS, "GL_UNPACK_SKIP_PIXELS"),
        (GL_UNPACK_SKIP_ROWS, "GL_UNPACK_SKIP_ROWS"),
        (GL_UNPACK_SWAP_BYTES, "GL_UNPACK_SWAP_BYTES"),
        (GL_VERTEX_ARRAY, "GL_VERTEX_ARRAY"),
        (GL_VERTEX_ARRAY_SIZE, "GL_VERTEX_ARRAY_SIZE"),
        (GL_VERTEX_ARRAY_STRIDE, "GL_VERTEX_ARRAY_STRIDE"),
        (GL_VERTEX_ARRAY_TYPE, "GL_VERTEX_ARRAY_TYPE"),
        (GL_VIEWPORT, "GL_VIEWPORT"),
        (GL_ZOOM_X, "GL_ZOOM_X"),
        (GL_ZOOM_Y, "GL_ZOOM_Y")]

    print
    for param, name in parameters:
        print "%-40s %s" % (name, glGetDoublev(param))
