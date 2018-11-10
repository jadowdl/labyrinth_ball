# This is a gigantic hack for now.
import sys
import math

LAT_DIVISIONS = 100
LNG_DIVISIONS = 100

# In a given fixation, all of the pathpoints should be within 
# DOT_DEGREES of each other on the unit sphere, great-circle-wise.
DOT_DEGREES = 58

# should be >= Level.java's HORIZON_LOCK_DOT_POINT
DOTLOCK = math.sin(math.pi * DOT_DEGREES/180.0)

NUM_CHIME_POINTS = 5

def dot(a, b):
    assert len(a) >= 3
    assert len(b) >= 3
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def distance_between(a,b):
    assert len(a) >= 3
    assert len(b) >= 3
    return math.sqrt(sum([(a[i]-b[i])**2 for i in range(0, 3)]))

def normalize(pt):
    mag = pt[0]*pt[0] + pt[1]*pt[1] + pt[2]*pt[2]
    assert(mag > 0.0)
    pt[0] /= mag
    pt[1] /= mag
    pt[2] /= mag

def firstZeroAfterPathDistance(points, target_distance):
    # get to distance
    i = 1
    d = 0.0
    while(i < len(points)):
        d += distance_between(points[i-1], points[i])
        if(d > target_distance and points[i][3] == 0):
            return i
        i += 1
    assert False

# on unit sphere
def get_point(lat_lin_coeff, lng_lin_coeff):
    assert math.fabs(lat_lin_coeff) <= 1.0
    assert math.fabs(lng_lin_coeff) <= 1.0

    lat_theta = math.pi * lat_lin_coeff
    lng_theta = -math.pi + 2.0*math.pi * lng_lin_coeff
    x = math.sin(lat_theta) * math.cos(lng_theta)
    y = math.sin(lat_theta) * math.sin(lng_theta)
    z = math.cos(lat_theta)    
    return [x,y,z]

# Latitude is represented as pi * t, for some linear coefficient t.  By setting
# minimum and maximum values of t, you can essentially cut off points near the
# poles as being candidates for camera focuses.
#
# Current values remove from consideration latitude less than pi/6 and greater
# than 5/6 * pi, which seems to do pretty well for our models.
MIN_LAT_COEFF = .16666
MAX_LAT_COEFF = .83333
assert MAX_LAT_COEFF > MIN_LAT_COEFF
def full_spel_set():
    ret = []
    for i in xrange(0, LAT_DIVISIONS):
        for j in xrange(0, LNG_DIVISIONS):
            true_lat_coeff = float(i)/LAT_DIVISIONS
            lat_coeff = MIN_LAT_COEFF + true_lat_coeff * (MAX_LAT_COEFF-MIN_LAT_COEFF)
            lng_coeff = float(j)/LNG_DIVISIONS
            ret.append(get_point( lat_coeff, lng_coeff)) 
    return ret

# Takes advantage of fact that all points in spels and point itself have unit length
def filter_spel_set(spels, point):
    ret = []
    for p in spels:
        if (dot(p, point) > DOTLOCK):
            ret.append(p)
    return ret

def printpoint(point):
    print "%.6f %.6f %.6f %d"%tuple(point)

def run():
    # parse input
    print "Parsing Input from file: %s"%sys.argv[1]
    points = []
    for l in open(sys.argv[1]).readlines():
        tokens = l.strip().split(" ")
        assert len(tokens) == 4
        point = [float(t) for t in tokens]
        normalize(point)
        # zero out path points - we're resetting and won't accept prior work.
        point[3] = 0
        # ignore camera focuses from before.  Add all others to the path.
        if (int(tokens[3]) != 3):
            points.append(point)
    print "  (Found %s path points)"%len(points)
    camfocuses = []

    # find fixations; for each set a cam trigger (type==2) and a cam focus (==3)
    print "Finding Fixations"
    i = 0
    j = 0
    while(i < len(points) and j < len(points)):
        print "  (i=%s)"%i
        focus_candidates = full_spel_set()
        filtered_focus_candidates = focus_candidates 
        j = i
        while(j < len(points) and len(filtered_focus_candidates) > 0):
            focus_candidates = filtered_focus_candidates
            filtered_focus_candidates = filter_spel_set(focus_candidates, points[j])
            j += 1
        # we have defined a fixation.  Set the cam trigger and focus for it.
        # we use focus_candidates[0], but we could use anything in the set.
        assert len(focus_candidates) > 0
        assert (j > i+1)
        # points[j-1][3] = 1  # hack - chime point, for debugging.
        camfocus = [k for k in focus_candidates[0]]
        camfocus.append(3) # type cam focus
        camfocuses.append(camfocus)

        # this fixations cam trigger becomes start point for next fixation
        i = i + int((j-i)/2)
        assert points[i][3] == 0
        points[i][3] = 2 # type cam trigger
    
    # Set Chime Points and move some triggers around.
    # 1) Set the first and last chime point
    assert points[0][3] == 0
    points[0][3] = 1 # type chime point
    assert points[-1][3] == 0
    points[-1][3] = 1 # type chim point
    print "Added goal chimes."

    # 2) find total path distance (cheat, just use Euclidean)
    path_distance = sum( map( lambda x: distance_between(x[0], x[1]), zip(points, points[1:])))
    print "Computed Total Path Distance As: ", path_distance

    # 3) Set NUM_CHIME_POINTS chime points in the middle, one every
    # approximately $path_distance / (NUM_CHIME_POINTS+1) along.
    for i in xrange(0, NUM_CHIME_POINTS):
        index = firstZeroAfterPathDistance(points, path_distance/(NUM_CHIME_POINTS+1)*(i+1))
        assert points[index][3] == 0
        print "Adding chimepoint at index: ", index
        points[index][3] = 1  # chime point

    # 4) Add an initial cam trigger
    DISTANCE_POLE_TO_TRIGGER = .25
    index = firstZeroAfterPathDistance(points, DISTANCE_POLE_TO_TRIGGER)
    print "Adding first cam trigger at index", index
    assert points[index][3] == 0
    points[index][3] = 2  # type cam trigger

    # 5) move the last cam trigger closer to the pole.
    for i in xrange(1, len(points)):
        if(points[-i][3] == 2):
            print "Deleting cam trigger at index:", -i
            points[-i][3] = 0
            break
    index = firstZeroAfterPathDistance(points, path_distance - DISTANCE_POLE_TO_TRIGGER)
    assert points[index][3] == 0
    print "Adding a new cam trigger at index:", index
    points[index][3] = 2  # type cam trigger
   
    # output
    print "==== NEW POINTS ===="
    for k in points:
        printpoint(k)
    for k in camfocuses:
        printpoint(k)

run()
