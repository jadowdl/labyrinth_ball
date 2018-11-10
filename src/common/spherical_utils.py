#!/usr/bin/env python2

#########################################
# (c) 2011 Evan Mallory, TSM Games
# (c) 2018 James Dowdell, TSM Games
#########################################

from math import sqrt, sin, pi, atan2, acos, cos, fabs
import numpy
import random # for unit test

# Takes a list [x, y, z] and returns [r, theta, phi]
# theta, the rotation around the x-axis, is in the range [0, pi].
# phi, the rotation around the y-axis, is in the range [-pi, pi].
# theta == "latitude"
# phi == "longitude"
def xyzToSphericalCoords(point):
    r = sqrt(sum([x*x for x in point]))
    theta = acos(point[2] / r)
    phi = atan2(point[1], point[0])
    return (r, theta, phi)


def sphericalCoordsToXyz(r, theta, phi):
    return (r * sin(theta) * cos(phi),
        r * sin(theta) * sin(phi),
        r * cos(theta))


def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step


def latLongToSinusoidal2d(latitude, longitude, width, height):
    x = int(.5 + (longitude * cos(latitude) * ((width - 1)/ (2 * pi)) + (width / 2.0)))
    y = int(.5 + ((latitude + pi/2) * ((height - 1) / pi)))
    if(x < 0): x = 0
    if(x >= width): x = width-1
    if(y < 0): y = 0
    if(y >= height): y = height-1
    return (x, y)


def sinusoidal2dToLatLong(x, y, width, height):
    latitude = y / ((height - 1) / pi) - pi/2
    longitude = (x - width/2.0)/((width - 1)/ (2 * pi))/cos(latitude)
    return (latitude, longitude)


def testOneSinusBckFwd(x, y, width, height, output=False):
    (lat, lng) = sinusoidal2dToLatLong(x, y, 800, 400)
    lat += pi / 2
    spherical = [sin(lat) * cos(lng),
                 sin(lat) * sin(lng),
                 cos(lat)]
    (newx, newy) = xyzToSinusoidal2d(spherical, 800, 400)

    if output:
        print "%d %d ====> %d %d" % (x, y, newx, newy)
    return (newx, newy)


def testSinusBckFwd():
    num_success = 0
    num_failure = 0
    num_x_failure = 0
    num_y_failure = 0
    num_both_failure = 0
    for x in range(150, 650):
        for y in range(100, 300):
            (newx, newy) = testOneSinusBckFwd(x, y, 800, 400)
            if x == newx and y == newy:
                num_success += 1
            else:
                num_failure += 1
                if x != newx and y == newy:
                    num_x_failure += 1
                if y != newy and x == newx:
                    num_y_failure += 1
                if x != newx and y != newy:
                    if num_both_failure < 50:
                        print "%d %d ==> %d %d" % (x, y, newx, newy)
                    num_both_failure += 1

    print "%d successes, %d failures" % (num_success, num_failure)
    print "%d x, %d y, %d both" % (num_x_failure, num_y_failure, num_both_failure)


def testSinusFwdBck():
    (lat, lng) = (pi/2 - pi*random.random(), pi - 2*pi*random.random())
    (x,y) = latLongToSinusoidal2d(lat, lng, 400, 200)
    (new_lat, new_lng) = sinusoidal2dToLatLong(x, y, 400, 200)
    #print lat, lng, x, y, new_lat, new_lng
    return (new_lat-lat, new_lng-lng)


# This actually breaks because of singularities at the poles, but it
# worked well enough that I'm satisfied that it works in general.
def unitTestSinus():
    for i in range(0, 100):
        t = testSinusFwdBck()
        assert fabs(t[0]) < .2 and \
               fabs(t[1]) < .2, "Margin of error too big: (%s, %s)"%(t)


def xyzToSinusoidal2d(p, width, height):
    (r, theta, phi) = xyzToSphericalCoords(p)
    latitude = theta - pi/2
    longitude = phi

    return latLongToSinusoidal2d(latitude, longitude, width, height)


# In "regular" land, (r, theta) in 2D <=> (pi/2 - pi*r/(h/2), theta) spherical
def latLongToRegular2d(lat, long, width, height):
    assert width == height
    assert -pi/2.0 <= lat and lat <= pi/2.0
    assert -pi <= long and long <= pi
    theta = long
    r = (height/2.0 - (lat/pi)*height) / 2.0
    x = r * cos(theta)
    x += width /2.0
    x *= float(width-1)/width # force it to be in range [0, width)
    y = r * sin(theta)
    y += height/2.0
    y *= float(height-1)/height # force it to be in range [0, height)
    return (x, y)


def regularRThetaToLatLong(r, theta, width, height):
    assert width == height
    assert r <= width*.5
    r /= (width*.5)
    long = theta
    lat = pi*.5 - pi*r
    return (lat, long)


def normalizePoint(p):
    a = numpy.array(p, 'f')
    return a / sqrt(numpy.dot(a, a))


if __name__ == "__main__":
    testSinusBckFwd()
