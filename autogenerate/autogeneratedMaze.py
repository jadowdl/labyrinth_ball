import sys
from PIL import Image
import math

sys.path.insert(0, '/Users/jdowdell/Downloads/python')
from makeModelFromSineMap import Point, loadPrecomputedPoints, displayPoints
from utils import regularRThetaToLatLong

###############################################################################
## CONSTANTS
###############################################################################

# circuit 0 == center / North Pole
# circuit 1 == first ring from center
# ...
# circuit NUM_CIRCUITS-1 == outtermost visible ring
# circuit NUM_CIRCUITS == outside / South Pole.
# Chartes Labyrinths has NUM_CIRCUITS = 12
NUM_CIRCUITS = 24

# For any (r, theta) points, we work in terms of the unit disc (so r is in 
# range [0,1]) and theta is in PRADIANS = Radians / (2pi)
# The below are all in PRADIANS = Radians / (2pi)
NORTH_POLE_WIDTH = .075# should be > 2*CIRCUIT_WIDTH, essentially size(circuit 0)
SOUTH_POLE_WIDTH = .075# should be > 0.  I suspect, for sphere, that it should
                      # also obey 2*CIRCUIT_WIDTH (or pseudo radials will be
                      # wider than it can handle that's what she said ohhh).
CIRCUIT_WIDTH = (1.0 - NORTH_POLE_WIDTH - SOUTH_POLE_WIDTH)/(NUM_CIRCUITS)
assert NORTH_POLE_WIDTH > 2*CIRCUIT_WIDTH
FIELD_WIDTH = .5 * CIRCUIT_WIDTH     
PATH_WIDTH = CIRCUIT_WIDTH - FIELD_WIDTH

IMAGE_SIZE = 1024

DISTORT_FOR_SPHERE = True

TYPE_NONE = 0
TYPE_LEFT = 1
TYPE_CENTER = 2
TYPE_RIGHT = 3
###############################################################################
## Helper Functions
###############################################################################
def CHECK_LIN_COEFF(t):
    assert t >= 0 and t <= 1.0

def CHECK_THREE_EQUAL(a,b,c):
    assert a==b and b == c

def MOD_ONE_POINT_OH(t, force_positive=False):
    i = int(t)
    t -= i
    if(force_positive and t < 0):
        t += 1.0
    return t

# point is (circuit_number, theta, k), where k is a float and represents
# signed distance from theta (positive goes towards the ccw side).  Note that
# we _translate_ by delta, not _rotate_, even though we're working in polar.
def confirm_point(pt):
    assert len(pt) == 3
    CHECK_LIN_COEFF(pt[0] / (NUM_CIRCUITS))
    assert(pt[0] == int(pt[0]))
    CHECK_LIN_COEFF(math.fabs(pt[1]))

# reverse the sign of (point1[1] - point2[1]) while keeping them in (-1, 1)
# requires the two thetas are < 1.0 apart (avoid nasty edge cases for  <=, like 
# start==0.0 and stop == 1.0
def reverse_orientation(point1, point2):
    confirm_point(point1)
    confirm_point(point2)
    assert math.fabs(point1[1]-point2[1]) < 1.0
    new_point1 = [i for i in point1]
    new_point1[1] += -1.0 * sign_of(point1[1]-point2[1])
    new_point2 = [i for i in point2]
    while(new_point1[1] <= -1.0):
        new_point1[1] += 1.0
        new_point2[1] += 1.0
    while(new_point1[1] >= 1.0):
        new_point1[1] -= 1.0
        new_point2[1] -= 1.0
    return (new_point1, new_point2)

# supports fraction circuit num; X+.5 is the center of circuit X
def CIRCUIT_NUM_TO_RADIUS(n):
    assert n >= 0.0 and n <= NUM_CIRCUITS
    return NORTH_POLE_WIDTH + n * CIRCUIT_WIDTH

# if DISTORT_FOR_SPHERE is false, this just adds theta_delta to theta and
# normalizes.  Otherwise, we 'peel' theta_delta away from theta according
# to the distortion formula (also found in mazeToMaze.py)
def sphere_adjustment(r, theta, theta_delta):
    CHECK_LIN_COEFF(r)
    if(DISTORT_FOR_SPHERE):
        theta_delta *= r * math.pi / \
            math.cos(math.pi/2 - r*math.pi)
    ret = MOD_ONE_POINT_OH(theta+theta_delta)
    return ret

# A ray "B" shoots out from the origin along theta.  Standing at the origin
# and looking down B, another ray "R" sits to its right.  R runs parallel
# to B, and at all times lies delta away from B.  We call R a pseudo radial.
# 
# R intersects the ring designated by ring_num in at most one point; iff delta
# is less than the ring's radius it won't intersect.  We assert on an
# intersection existing, and return the intersection point in (r, theta)
#
# The deriviation of theta is easy if you rotate the circle so that B lies at 0
# degrees.  You rise or fail delta from the x-axis, and any trig class gives
# the answer - its just arcsin(delta) away from B.
def intersectPseudoRadialWithRing(theta, delta, ring_radius):
    assert math.fabs(delta) <= ring_radius
    CHECK_LIN_COEFF(ring_radius)  # implicitly checks delta as well, via above
    CHECK_LIN_COEFF(math.fabs(theta))
    theta_delta = math.asin(delta / ring_radius) / (2*math.pi) # pi/2 / 4 = pradians!
    return (ring_radius, sphere_adjustment(ring_radius, theta, theta_delta))

###############################################################################
## CLASSES
###############################################################################

# poor man's abstract base class for our three guide types.
class Guide:
    # start and stop have same format: (circuit_num, theta)
    def getPathLCR(self, start, stop, t):
        assert False
    
    def getPathLCRWithParity(self, start, stop, t, parity=False):
        r = self.getPathLCR(start, stop, t)
        if(parity):
            return r
        else:
            return (r[2], r[1], r[0])

class RingGuide(Guide):
    def __init__(self, circuit_num):
        CHECK_LIN_COEFF(circuit_num / float(NUM_CIRCUITS))
        self.circuit_num = circuit_num

    def getPathLCR(self, start, stop, t):
        CHECK_LIN_COEFF(t)
        CHECK_THREE_EQUAL(self.circuit_num, start[0], stop[0])
        CHECK_LIN_COEFF(math.fabs(start[1]))
        CHECK_LIN_COEFF(math.fabs(stop[1]))
        r = CIRCUIT_NUM_TO_RADIUS(.5 + self.circuit_num)
        theta = start[1] + t * (stop[1] - start[1])
        delta = .5 * PATH_WIDTH
        return ((r-delta, theta), (r, theta), (r+delta, theta))

# This class assumes all start and stop points lie on cardinal lines
# (theta = 0, .25, .5, .75, or 1.0).  To prevent overshooting a sweep
# along the normal RingGuide and line up perfectly with a PivotCircle 
# lying near the cardinal, margins are introduced by adjusting the
# actual theta passed on to RingGuide.getPathLCR.
class AutoAdjustingRingGuide(RingGuide):
    def __init__(self, circuit_num, extra_southern_offset=True):
        RingGuide.__init__(self, circuit_num)
        self.extra_southern_offset = extra_southern_offset

    # Ideally, we would perform adjustments for L, C, and R; for now
    # just adjust for C, and assume it's "close enough" for L and R.
    def _adjust(self, point, cw_orientation=True):
        ring_radius = CIRCUIT_NUM_TO_RADIUS(point[0]+.5)
        delta = CIRCUIT_WIDTH
        # southern cardinal gets extra margin because everything's pushed
        # outwards by the visible pseudo-radial lines.
        if(self.extra_southern_offset and \
            (point[1] == .75 or point[1] == -.25)):
            delta += CIRCUIT_WIDTH
        if(not cw_orientation):
            delta *= -1.0
        (dummy, new_theta) =  \
            intersectPseudoRadialWithRing(point[1], delta, ring_radius)
        return (point[0], new_theta)

    def getPathLCR(self, start, stop, t):
        assert math.fabs(start[1]) in set([0.0, .25, .5, .75, 1.0])
        assert math.fabs( stop[1]) in set([0.0, .25, .5, .75, 1.0])
        return RingGuide.getPathLCR(self, 
                                self._adjust(start, start[1] < stop[1]), 
                                self._adjust(stop,  start[1] > stop[1]), 
                                t)

# Experimental - less automatic than AutoAdjustingRingGuide
# Also enforces that the direction of travel agrees with start's delta
class AdjustingRingGuide(RingGuide):
    # Ideally, we would perform adjustments for L, C, and R; for now
    # just adjust for C, and assume it's "close enough" for L and R.
    def _adjust(self, point, delta):
        ring_radius = CIRCUIT_NUM_TO_RADIUS(point[0]+.5)
        (dummy, new_theta) =  \
            intersectPseudoRadialWithRing(point[1], delta, ring_radius)
        return (point[0], new_theta)

    # TODO(jdowdell) - this is weak against start[2] == 0
    # in that case, consider just checking that we dont sweep > .5
    def getPathLCR(self, start, stop, t):
        if(sign_of(stop[1]-start[1]) == -sign_of(start[2])):
            (start, stop) = reverse_orientation(start, stop)
        return RingGuide.getPathLCR(self, 
                                self._adjust(start, start[2]*CIRCUIT_WIDTH), 
                                self._adjust(stop,   stop[2]*CIRCUIT_WIDTH), 
                                t)

class PseudoRadialGuide(Guide):
    # theta gives the True Radial
    # delta gives the offset to make the pseudo radial.  See comments for
    # intersectPseudoRadialWithRing()
    def __init__(self, theta, delta):
        CHECK_LIN_COEFF(theta)
        CHECK_LIN_COEFF(math.fabs(delta))
        self.theta = theta
        self.delta = delta

    # [R]adius [To] [I]ntersection (Point)
    def r2I(self, radius, margin):
        return intersectPseudoRadialWithRing(self.theta, 
            self.delta + margin, radius)

    def getPathLCR(self, start, stop, t):
        CHECK_LIN_COEFF(t)
        CHECK_THREE_EQUAL(self.theta, start[1], stop[1])
        CHECK_LIN_COEFF(start[0]/NUM_CIRCUITS)
        CHECK_LIN_COEFF(stop[0]/NUM_CIRCUITS)
        # This takes advantage of the fact that all pseudo radials we use run
        # origin-ward.  The innermost stretch of stop's circuit is given by stop[0],
        # the outtermost stretch of start's circuit is actually the innermost of
        # (start[0]+1)'s
        start_r = CIRCUIT_NUM_TO_RADIUS(1 + start[0])
        stop_r  = CIRCUIT_NUM_TO_RADIUS(stop[0])
        current_r = start_r + t * (stop_r - start_r)
        margin = .5*PATH_WIDTH
        return (self.r2I(current_r, -margin),
                self.r2I(current_r, 0.0),
                self.r2I(current_r, margin))

class PivotCircleGuide(Guide):
    def __init__(self, circuit_num, true_radial_theta, pseudo_radial_delta):
        #TODO
        # print "DEBUG", circuit_num,",",true_radial_theta,",",pseudo_radial_delta
        CHECK_LIN_COEFF(circuit_num / float(NUM_CIRCUITS))
        CHECK_LIN_COEFF(true_radial_theta)
        self.circuit_num = circuit_num
        self.true_radial_theta = true_radial_theta
        # pivot sits in the field at the juncture between two circuits.
        # it will join circuit_num with (circuit_num-1)
        circuit_radius = CIRCUIT_NUM_TO_RADIUS(self.circuit_num)
        (self.center_r, self.center_theta) = intersectPseudoRadialWithRing(\
            true_radial_theta, pseudo_radial_delta, circuit_radius)

    # [m]inor [To] [M]ajor, ie convert off origin polar circle (minor) treating 
    # self.center_r and self.center_theta as the origin, to regular (major) 
    # polar coords.
    #
    # I'm currently cheating and using an extremely naive and expensive 
    # version, because middle school geometry was kicking my ass.
    # 
    # Of potential use in the future:
    #   1) Equation of off-origin circle: http://en.wikipedia.org/wiki/Circle
    #   2) I'm pretty sure that if you simplify the below, major_r becomes 
    #      math.sqrt(self.center_r**2 + CIRCUIT_WIDTH**2 + 
    #                   2 * self.center_r * CIRCUIT_WIDTH * math.cos(
    #                   2.0 * math.pi * (self.center_theta - minor_theta) ))
    #   3) given major_r, you have classic a Side/Side/Side triangle setup, 
    #   and should be able to add to self.center_theta accordingly.
    def m2M(self, minor_r, minor_theta):  
        x = self.center_r * math.cos(2.0*math.pi*self.center_theta) +\
            minor_r * math.cos(2.0*math.pi*minor_theta)
        y = self.center_r * math.sin(2.0*math.pi*self.center_theta) +\
            minor_r * math.sin(2.0*math.pi*minor_theta)
        major_r = math.sqrt(x**2 + y**2)
        major_theta = MOD_ONE_POINT_OH(math.atan2(y,x)/(2.0*math.pi))
        # if it weren't for sphere adjustments, we'd just return here.
        delt = MOD_ONE_POINT_OH(major_theta - self.center_theta)
        if (math.fabs(delt) > .5): # always take shortest distance
            delt -= sign_of(delt)
        assert(math.fabs(delt) < 1.0), "%s %s %s"%(major_theta, self.center_theta, delt)
        return (major_r, sphere_adjustment(major_r, self.center_theta, delt))
    
    # Ultimately very similar to RingGuide
    def getPathLCR(self, start, stop, t):
        CHECK_LIN_COEFF(t)
        CHECK_THREE_EQUAL(self.circuit_num, start[0], stop[0])
        CHECK_LIN_COEFF(math.fabs(start[1]))
        CHECK_LIN_COEFF(math.fabs(stop[1]))
        minor_theta = start[1] + t * (stop[1] - start[1])
        center = CIRCUIT_WIDTH/2
        inside_margin = PATH_WIDTH/2
        return (self.m2M(center-inside_margin, minor_theta),
                self.m2M(center,               minor_theta),
                self.m2M(center+inside_margin, minor_theta))

def sign_of(x):
    if (x<0):
        return -1.0
    elif (x > 0):
        return 1.0
    return 0.0

# PivotCircleGuide expects the thetas in start and stop are relative to it's center,
# not the origin.  This class translates to what PivotCircleGuide expects.
class LiteralPivotCircleGuide(PivotCircleGuide):
    def getPathLCR(self, start, stop, t):
        # TODO - verify start and stop
        confirm_point(start)
        confirm_point(stop)
        assert(stop[2]-start[2] in set([-.5, 0, .5]))
        assert(sign_of(start[2]) == sign_of(stop[2]))
        assert(math.fabs(start[0] - stop[0]) == 1)
        # assert that start[1] equivalent to stop[1]?

        # sweep distance is in pradians; semicircle if deltas agree, generally 90 degrees otherwise.
        # if positive we sweep ccw from start; if negative, cw
        sweep_distance = .5 * sign_of(start[2]) * sign_of(stop[0]-start[0])

        # We're basically going to create a new start (start_copy) and new
        # stop (stop_copy) from the ground up.  Start with 'start' though.
        start_copy = [i for i in start]
        # Offset theta bit so that entry and exit into the pivot circle agrees
        # with the ring guide abutting it.
        start_copy[1] += self.center_theta - self.true_radial_theta
        # Make sure we have the right circuit...
        if(start_copy[0] != self.circuit_num):
            assert(stop[0] == self.circuit_num)
            start_copy[0] = stop[0]
        # Possibly Reverse start_copy's theta...
        if (start[0] < stop[0]): # inner circuit lies .5 pradians away from radial orientation.
            start_copy[1] = MOD_ONE_POINT_OH(start_copy[1] + .5)
        # Base stop_copy off start_copy, but add sweep_distance to the theta
        stop_copy = [i for i in start_copy]
        stop_copy[1] = stop_copy[1]+sweep_distance
        # treat misaligned deltas as though halfway through a normal pivot
        if(start[2] != stop[2]):
            sweep_distance /= 2;
            if(math.fabs(start[2]) < math.fabs(stop[2])):
                start_copy[1] = start_copy[1] + sweep_distance
            else:
                stop_copy[1] = stop_copy[1] - sweep_distance
        # MOD_ONE_POINT_OH can distort orientation of start -> stop; so we
        # have to manually enforce thetas in (-1, 1).  We know we started with
        # start in (-1, 1), and math.fabs(sweep_distance) < .5, so both of these
        # should be in (-2, 2), and no more than .5 apart.
        while(start_copy[1] <= -1.0 or stop_copy[1] <= -1.0):
            start_copy[1] += 1.0
            stop_copy[1] += 1.0
        while(start_copy[1] >= 1.0 or stop_copy[1] >= 1.0):
            start_copy[1] -= 1.0
            stop_copy[1] -= 1.0
        CHECK_LIN_COEFF(math.fabs(start_copy[1]))
        CHECK_LIN_COEFF(math.fabs(stop_copy[1]))

        # print "DEBUG Translated %s => %s, %s => %s"%(start, start_copy, stop, stop_copy)
        return PivotCircleGuide.getPathLCR(self, start_copy, stop_copy, t)

###############################################################################
## ENTRY CODE 
###############################################################################

# (r, theta) is on unit disc, theta is in pradians
# basically mapping to a square of side image_width
# We flip y, assuming that we're about to be used on an image.
def rThetaToXY( pt , image_width, toInt = False):
    (r, theta) = pt
    assert -1.0 <= r
    assert  r <= 1.0 
    actual_radius = image_width/2.0
    x = image_width/2.0
    y = image_width/2.0
    x += actual_radius * r * math.cos(theta * 2.0 * math.pi)
    y -= actual_radius * r * math.sin(theta * 2.0 * math.pi)
    x *= float(image_width-1)/image_width
    y *= float(image_width-1)/image_width
    if(toInt):
        x = int(round(x))
        y = int(round(y))
    return (x,y)

def apply_callback_on_path_seg(callback, guide, start, stop, parity=False, STEPS=1000):
   for i in xrange(0, STEPS):
       pts = guide.getPathLCRWithParity(start, stop, i/float(STEPS), parity)
       apply(callback, (pts[0], TYPE_LEFT))
       apply(callback, (pts[1], TYPE_CENTER))
       apply(callback, (pts[2], TYPE_RIGHT))

def draw_on_image_callback(image, pt, type):
    assert type in set([TYPE_LEFT, TYPE_CENTER, TYPE_RIGHT])
    color = (255, 255, 255)
    if (type == TYPE_LEFT):
        color = (255, 255, 0)
    elif (type == TYPE_RIGHT):
        color = (0, 255, 0)
    image.putpixel(rThetaToXY(pt, IMAGE_SIZE, True), color)
    
def debug_circuits(image):
    STEPS = 10000
    for k in range(0, NUM_CIRCUITS):
        g = RingGuide(k)
        start = (k, 0.0)
        stop  = (k, 1.0)
        draw_on_image(image, g, start, stop)

def debug_radials(image):
    STEPS = 1000
    SWEEP = 4
    for k in range(0, SWEEP):
        theta = k * 1.0/SWEEP
        delta = CIRCUIT_WIDTH
        g = PseudoRadialGuide(theta, delta)
        start = (1              , theta)
        stop  = (NUM_CIRCUITS-1 , theta)
        draw_on_image(image, g, start, stop)

def debug_pivots(image):
    STEPS = 1000
    pivot_specs = [0  , -  CIRCUIT_WIDTH, # NORTH OF EAST
                   .25,    CIRCUIT_WIDTH, # EAST OF NORTH
                   .25, -  CIRCUIT_WIDTH, # WEST OF NORTH
                   .5 , -  CIRCUIT_WIDTH, # NORTH OF WEST
                   .5 ,    CIRCUIT_WIDTH, # SOUTH OF WEST
                   .75, -2*CIRCUIT_WIDTH, # WEST OF SOUTH
                   .75, -  CIRCUIT_WIDTH, # VISIBLE PSEUDO WEST OF SOUTH
                   .75,    CIRCUIT_WIDTH, # VISIBLE PSEUDO EAST OF SOUTH
                   .75,  2*CIRCUIT_WIDTH, # EAST OF SOUTH
                   0  ,   CIRCUIT_WIDTH]  # SOUTH OF EAST
    for a in range(0, 10):
        for k in range(2, NUM_CIRCUITS-1): # goes off image at outtermost.
            g = PivotCircleGuide(k, pivot_specs[a*2], pivot_specs[1+2*a])
            start = (k, 0.0)
            stop = (k, 1.0)
            draw_on_image(image, g, start, stop)


def manual_path(output_image):
    # I can condense this to just an array of (n, theta) points, and after
    # that generalize the Chartes formula so manual specification is
    # unnecessary.
    #
    # KNOWN BUGS (letters) / POINTS OF CONFUSION (lower case romans)
    #  i) Visible Pseudo Radial lines are offset by CIRCUIT_WIDTH/2
    #  though in many of my notes I had them listed as offset by
    #  a full CIRCUIT_WIDTH
    #  ii) PivotCircleGuide circuit numbers are sometimes off by 1 from
    #  what you would expect
    #  A) All guides overshoot because they are interpreting me literally; 
    #  they actually need to cut in later and out sooner than the cardinal
    #  of 0, .25, .5, or .75
    #  B/iii) Yellow and Green naturally switch back and forth, 
    #  even though they should be connected.  We need to incorporate some
    #  concept of handed-ness or orientation.  Currently manually fixed via
    #  specifying parity to draw_on_image; i can't see the pattern either.
    #  iv) Deltas for pseudo radials are not absolute, but twist as we go 
    #  around.  I should not be specifying the same sign
    #  on the delta for both south and north radials; rather the delta
    #  should always go 'CCW' or 'CW', and so inverts sign in that case.
    #  C) Theta specifications should be allowed to be in -1 to 1, instead
    #  of 0 to 1, to allow for going CW vs. CCW along a guide.  To swing
    #  through the x-axis, currently have to draw twice (!)
    g = PseudoRadialGuide(.75, -CIRCUIT_WIDTH/2)
    draw_on_image(output_image, g, (11, 0.75), (7, 0.75), False)
    g = PivotCircleGuide(8, 0.75, -CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (8, 0.0), (8,0.25), False)
    g = AutoAdjustingRingGuide(7, False) # False bc/ exiting pseudo-radial
    draw_on_image(output_image, g, (7, 0.75), (7, 0.5), True)
    g = PivotCircleGuide(7, 0.50, +CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (7, 0.5), (7, 0.0), True)
    g = AutoAdjustingRingGuide(6, False) # False bc/ entering pseudo-radial 
    draw_on_image(output_image, g, (6, 0.50), (6, .75), False)
    g = PivotCircleGuide(6, 0.75, -CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (6, .75), (6, 1.0), False)
    g = PseudoRadialGuide(.75, -CIRCUIT_WIDTH/2)
    draw_on_image(output_image, g, (6, 0.75), (1, 0.75), False)
    g = PivotCircleGuide(2, 0.75, -CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (2, .0), (2, .25), False)
    g = AutoAdjustingRingGuide(1, False) # False bc/ exiting pseudo-radial
    draw_on_image(output_image, g, (1, 0.75), (1, 0.25), True)
    g = PivotCircleGuide(2, 0.25, +CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (2, -.25), (2, .25), False)
    g = AutoAdjustingRingGuide(2)
    draw_on_image(output_image, g, (2, 0.25), (2, 0.75), False)
    g = PivotCircleGuide(3, 0.75, -2*CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (3, .25), (3, -.25), True)
    g = AutoAdjustingRingGuide(3)
    draw_on_image(output_image, g, (3, 0.75), (3, 0.50), True)
    g = PivotCircleGuide(4, 0.50,  CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (4, .0), (4, .5), False)
    g = AutoAdjustingRingGuide(4)
    draw_on_image(output_image, g, (4, 0.5), (4, 0.75), False)
    g = PivotCircleGuide(5, 0.75,  -2*CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (5, .25), (5, -.25), True)
    g = AutoAdjustingRingGuide(5)
    draw_on_image(output_image, g, (5, 0.75), (5, 0.25), True)
    g = PivotCircleGuide(5, 0.25,  CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (5, .25), (5, -.25), True)
    g = AutoAdjustingRingGuide(4)
    draw_on_image(output_image, g, (4, 0.25), (4, 0.5), False)
    g = PivotCircleGuide(4, 0.50,  -CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (4, .5), (4, 1.0), False)
    g = AutoAdjustingRingGuide(3)
    draw_on_image(output_image, g, (3, 0.50), (3, 0.00), True)
    g = PivotCircleGuide(3, 0.00,  +CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (3, 1.0), (3, 0.50), True)
    g = AutoAdjustingRingGuide(2)
    draw_on_image(output_image, g, (2, 0.00), (2, 0.25), False)
    g = PivotCircleGuide(2, 0.25,  -CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (2, .25), (2, 0.75), False)
    g = AutoAdjustingRingGuide(1)
    draw_on_image(output_image, g, (1, 0.25), (1, -.25), True)
    g = PivotCircleGuide(2, 0.75,  2*CIRCUIT_WIDTH) # South == ++
    draw_on_image(output_image, g, (2, .25), (2, 0.75), False)
    g = AutoAdjustingRingGuide(2)
    draw_on_image(output_image, g, (2, -.25), (2, .00), False)
    g = PivotCircleGuide(3, 0.00,  -CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (3, .50), (3, 0.00), True)
    g = AutoAdjustingRingGuide(3)
    draw_on_image(output_image, g, (3, 1.00), (3, .75), True)
    g = PivotCircleGuide(4, 0.75,  2*CIRCUIT_WIDTH) # South == ++
    draw_on_image(output_image, g, (4, .25), (4, 0.75), False)
    g = AutoAdjustingRingGuide(4)
    draw_on_image(output_image, g, (4, -.25), (4, .25), False)
    g = PivotCircleGuide(5, 0.25,  -CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (5, .75), (5, 0.25), True)
    g = AutoAdjustingRingGuide(5)
    draw_on_image(output_image, g, (5, .25), (5, .00), True)
    g = PivotCircleGuide(6, 0.00,  +CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (6, .5), (6, 1.00), False)
    g = AutoAdjustingRingGuide(6)
    draw_on_image(output_image, g, (6, .00), (6, .50), False)
    g = PivotCircleGuide(7, 0.50,  -CIRCUIT_WIDTH)
    draw_on_image(output_image, g, (7, 1.00), (7, .5), True)
    g = AutoAdjustingRingGuide(7)
    draw_on_image(output_image, g, (7, .50), (7, .25), True)


def choose_guide(start, stop):
    confirm_point(start)
    confirm_point(stop)
    g = None
    # CASE: Ring Guide - any time you start and stop on the same circuit.
    if (start[0] == stop[0] and
        MOD_ONE_POINT_OH(start[1]-stop[1]) != 0.0):
        g = AdjustingRingGuide(start[0])
    # Case: Pseudo Radial Guide (PRG); same theta, same delta, circuits more
    # than 1 apart.
    elif (MOD_ONE_POINT_OH(start[1]) == MOD_ONE_POINT_OH(stop[1]) and \
        start[2] == stop[2] and \
        math.fabs(start[0]-stop[0]) > 1):
        g = PseudoRadialGuide(start[1], start[2]*CIRCUIT_WIDTH)
    # CASE: Pivot Circle.  PRG -> RING = 90 degree arc; RING -> RING = 180
    # degree arc.  Circuits one apart, same theta
    # CASE: PRG -> Ring connector; circuits one apart, same theta, deltas
    # .5 apart.
    elif (math.fabs(start[0] - stop[0]) == 1 and \
          MOD_ONE_POINT_OH(start[1]) == MOD_ONE_POINT_OH(stop[1]) and \
          math.fabs(start[2] - stop[2]) in set([0.0, .5])):
        delta = start[2] 
        if(math.fabs(start[2]) < math.fabs(stop[2])): # set delta to the largest by magnitude
            delta = stop[2]
            # pivot circle at n links n -> (n-1), so go with max of start[0] vs. stop[0]
        g = LiteralPivotCircleGuide(max(start[0], stop[0]), start[1], delta*CIRCUIT_WIDTH)
    # CASE: UNKNOWN!
    else:
        assert False
    print start, ",", stop, " => ", g, "(", max(start[0], stop[0]),")"
    return g

# Walks the path, calling callback for points it touches.
#
# 'callback' takes arguments ((r, theta), type), where type is one of TYPE_LEFT
# TYPE CENTER or TYPE RIGHT.
# 
# parity switching happens before pivot circles if we're traveling away from
# the origin on a turn, and after them if we're traveling towards the origin on
# the turn.  Exception of pseudo radials:
#   1) resets parity to an xor of which side you entered from (global ccw or global
# cw) and what parity was
# winding orientation to horizontal delta position.
#   2) pivot-radial-pivot should be seen as one big pivot, so we need to restore
# to the opposite of the old parity when we're done, so that the second pivot
# recovers the effect the first would have had if it magically linked with the next
# ring guide.
def apply_callback_on_path(callback, path):
    parity = False
    tmp_parity = True 
    for i in xrange(1, len(path)):
        cur  = path[i]
        last = path[i-1]
        g = choose_guide(last, cur)  # also confirms cur and last are valid.
        if(isinstance(g, PseudoRadialGuide)):
            tmp_parity = parity
            parity = (last[2] < 0) != (parity) 
        elif(cur[0]>last[0] and isinstance(g, PivotCircleGuide)):
            parity = not parity
        apply_callback_on_path_seg(callback, g, last, cur, parity)
        if(isinstance(g, PseudoRadialGuide)):
            parity = not tmp_parity
        elif(last[0]>cur[0] and isinstance(g, PivotCircleGuide)):
            parity = not parity

# The beauty of the labyrinth ball is that it is the same maze whether you
# start from the north pole or the south.  This is because the maze diagram
# exhibits radial symmetry.  'Inverting' the first half of the maze yields the
# second half; and inverting any path segment X away from the beginning yields 
# the path segment X from the end.  Thus, many of Jensen's pieces can be 
# inverted into each other.  Explicitly define one, and you get the other free.
#
# 'Inversion' is inverted in just about every way you can imagine:
#    1) invert sign of all deltas
#    2) flip all thetas across y axis
#    3) invert circuit order
#    4) When you're finally done, reverse the path :)
def invert(path_segment):
    path_segment = [ (NUM_CIRCUITS-i[0], 
                      MOD_ONE_POINT_OH(.5 - i[1], True), 
                      -i[2]) for i in path_segment]
    path_segment.reverse()
    return path_segment

def jensen_outter_entry():
    ret = []
    n = NUM_CIRCUITS-1
    while(n > 0):
        ret.append((n  , .75, - .5))
        ret.append((n-3, .75, - .5))
        ret.append((n-4, .75, -1.0))
        if(n > 5):
            ret.append((n-4, .50,  1.0))
            ret.append((n-5, .50,  1.0))
            ret.append((n-5, .75, -1.0))
        n -= 6
    
    return ret

# This segment is essentially common to Jensen's L, M, and N
# It has its own self symmetry, so could probably be condensed ever farther;
# But I'll leave that fun for he who comes after me.  (probably me again.)
def jensen_critical_segment(n):
    return [\
    ( n+1, .25,  1.0),
    ( n+1, .75, -2.0),
    ( n+2, .75, -2.0),
    ( n+2, .50,  1.0),
    ( n+3, .50,  1.0),
    ( n+3, .75, -2.0),
    ( n+4, .75, -2.0),
    ( n+4, .25,  1.0),
    ( n+3, .25,  1.0),
    ( n+3, .50, -1.0),
    ( n+2, .50, -1.0),
    ( n+2, .00,  1.0),
    ( n+1, .00,  1.0),
    ( n+1, .25, -1.0),
    ( n+0, .25, -1.0),
    ( n+0, .75,  2.0),
    ( n+1, .75,  2.0),
    ( n+1, .00, -1.0),
    ( n+2, .00, -1.0),
    ( n+2, .75,  2.0),
    ( n+3, .75,  2.0),
    ( n+3, .25, -1.0),
    ( n+4, .25, -1.0),
    ( n+4, .00,  1.0),
    ( n+5, .00,  1.0),
    ( n+5, .25,  0.0)]

def jensen_N():
    ret = [(1, .25,  1.0)]
    ret.extend(jensen_critical_segment(1))
    return ret

def jensen_L():
    return invert(jensen_N())

def jensen_M(n):
    ret = [(n  , .25,  0.00),
           (n  , .50, -1.00),
           (n+1, .50, -1.00),
           (n+1, .25,  1.00)]
    ret.extend(jensen_critical_segment(n+1))
    return ret

def jensen_inner_entry():
    return invert(jensen_outter_entry())

def jensen_algorithm():
    assert(NUM_CIRCUITS >= 12 and NUM_CIRCUITS%6==0)
    path = []
    path.extend(jensen_outter_entry())
    n = 1 
    path.extend(jensen_N())
    n += 5
    while(n < NUM_CIRCUITS - 6):
        path.extend(jensen_M(n))
        n += 6
    path.extend(jensen_L())
    path.extend(jensen_inner_entry())
    # Every junction between two L, N, or M's causes a duplicate path point 
    # that should be removed.
    return remove_duplicates(path)

def remove_duplicates(path):
    new_path = []
    last = (-1, 0.0, 0.0)
    for p in path:
        if(p[0] == last[0] and \
           MOD_ONE_POINT_OH(p[1], True) == MOD_ONE_POINT_OH(last[1], True) and\
           math.fabs(p[2] - last[2]) < .0001):
            continue
        new_path.append(p)
        last = p
    return new_path

# "semi-manual" cause I had to figure out and specify all these points by hand :)
def semi_manual_path():
    return [\
    (11, .75, - .5),
    ( 8, .75, - .5),
    ( 7, .75, -1.0),
    ( 7, .50,  1.0),
    ( 6, .50,  1.0),
    ( 6, .75, -1.0),
    ( 5, .75, - .5),
    ( 2, .75, - .5),
    ( 1, .75, -1.0),
    ( 1, .25,  1.0),
    ( 2, .25,  1.0),
    ( 2, .75, -2.0),
    ( 3, .75, -2.0),
    ( 3, .50,  1.0),
    ( 4, .50,  1.0),
    ( 4, .75, -2.0),
    ( 5, .75, -2.0),
    ( 5, .25,  1.0),
    ( 4, .25,  1.0),
    ( 4, .50, -1.0),
    ( 3, .50, -1.0),
    ( 3, .00,  1.0),
    ( 2, .00,  1.0),
    ( 2, .25, -1.0),
    ( 1, .25, -1.0),
    ( 1, .75,  2.0),
    ( 2, .75,  2.0),
    ( 2, .00, -1.0),
    ( 3, .00, -1.0),
    ( 3, .75,  2.0),
    ( 4, .75,  2.0),
    ( 4, .25, -1.0),
    ( 5, .25, -1.0),
    ( 5, .00,  1.0),
    ( 6, .00,  1.0),
    ( 6, .25,  0.0), # midpoint.
    ( 6, .50, -1.0),
    ( 7, .50, -1.0),
    ( 7, .25,  1.0),
    ( 8, .25,  1.0), 
    ( 8, .75, -2.0),
    ( 9, .75, -2.0),
    ( 9, .50,  1.0),
    (10, .50,  1.0),
    (10, .75, -2.0),
    (11, .75, -2.0),
    (11, .25,  1.0),
    (10, .25,  1.0),
    (10, .50, -1.0),
    ( 9, .50, -1.0),
    ( 9, .00,  1.0),
    ( 8, .00,  1.0),
    ( 8, .25, -1.0),
    ( 7, .25, -1.0),
    ( 7, .75,  2.0),
    ( 8, .75,  2.0),
    ( 8, .00, -1.0),
    ( 9, .00, -1.0),
    ( 9, .75,  2.0),
    (10, .75,  2.0), 
    (10, .25, -1.0),
    (11, .25, -1.0),
    (11, .75,  1.0),
    (10, .75,   .5),
    ( 7, .75,   .5),
    ( 6, .75,  1.0),
    ( 6, .00, -1.0),
    ( 5, .00, -1.0),
    ( 5, .75,  1.0),
    ( 4, .75,   .5),
    ( 1, .75,   .5)]

def generate():
    assert len(sys.argv) == 2,\
        "Usage: %s <output_image>" % sys.argv[0]

    # Get Ready
    output_image_file = sys.argv[1]

    # Setup output image
    output_image = Image.new("RGB", (IMAGE_SIZE, IMAGE_SIZE), (0, 0, 0))
    
    # DEBUG
    # debug_circuits(output_image)
    # debug_radials(output_image)
    # debug_pivots(output_image)

    # MANUAL PATH
    # manual_path(output_image)    

    # "AUTO" PATH
    #    ATTEMPT 1
    ## assert NUM_CIRCUITS == 12
    ## path = semi_manual_path()

    #    ATTEMPT 2
    path = jensen_algorithm()

    apply_callback_on_path(\
        lambda x,y: draw_on_image_callback(output_image, x, y), 
        path)

    # Finish
    output_image.save(output_image_file)

###############################################################################
# LINK WITH makeModelFromSineMap.py
###############################################################################

def getKindForType(type):
    if(type == TYPE_LEFT):
        return Point.KIND_GREEN
    elif(type == TYPE_RIGHT):
        return Point.KIND_YELLOW
    elif(type == TYPE_CENTER):
        return Point.KIND_ENDPOINT
    else:
        assert False

# create a point of the variety found in makeModelFromSineMap
def createEvanPoint(r, theta, type):
    image_point_2d = rThetaToXY((r, theta), IMAGE_SIZE, True)

    r = r * IMAGE_SIZE/2
    theta = theta * math.pi * 2
    (latitude, longitude) = regularRThetaToLatLong(r, theta, IMAGE_SIZE, IMAGE_SIZE)
    
    assert( latitude <= math.pi/2 and latitude >= -math.pi/2 and
            longitude <= math.pi and longitude >= -math.pi)
    latitude += pi/2
    
    point = numpy.array([sin(latitude) * cos(longitude),
                         sin(latitude) * sin(longitude),
                         cos(latitude)], 'f')
    
    return Point(point, getKindForType(type), image_point_2d)

def append_points(startPoints, leftPoints, rightPoints, point, type):
    arr = None
    if (type == TYPE_CENTER and len(startPoints) == 0):
        arr = startPoints
    elif (type == TYPE_LEFT):
        arr = leftPoints
    elif (type == TYPE_RIGHT):
        arr = rightPoints

    if (arr is not None):
        arr.append(createEvanPoint(point[0], point[1], type))

def point_link(a, b):
    a.visited = True
    a.neighbors.append(b)

# really, compute neighbors
def finalize_points(startPoints, leftPoints, rightPoints):
    assert len(startPoints) > 0
    startPoints = startPoints[:1]
    if (len(leftPoints) > 0):
        point_link(startPoints[0], leftPoints[0])
    for i in xrange(1, len(leftPoints)):
        point_link(leftPoints[i-1], leftPoints[i])
    if (len(rightPoints) > 0):
        point_link(startPoints[0], rightPoints[0])
    for i in xrange(1, len(rightPoints)):
        point_link(rightPoints[i-1], rightPoints[i])

def generateModel():
    path = jensen_algorithm()
    startPoints = []
    leftPoints = []
    rightPoints = []
    apply_callback_on_path(\
        lambda x,y: append_points(startPoints, leftPoints, rightPoints, x, y),
        path)
    finalize_points(startPoints, leftPoints, rightPoints)
    precomputedPoints = []
    precomputedPoints.extend(startPoints)
    precomputedPoints.extend(rightPoints)
    precomputedPoints.extend(leftPoints)

    loadPrecomputedPoints(precomputedPoints)
    displayPoints

###############################################################################
# MAIN
###############################################################################

generateModel()
