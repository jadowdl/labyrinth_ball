#########################################
# (c) 2011 James Dowdell, TSM Games
#

import Image
import sys
from utils import *
from math import pi, cos, fabs
import math

# x in range 0..width, y in 0..height (inclusive!)
def lookupAliasedPixel(x, y, image):
    (w,h) = image.size
    x *= float(w-1)/w
    y *= float(h-1)/h
    horCoeff = x - int(x)
    verCoeff = y - int(y)
    base = [int(x), int(y)]
    dx = 1
    if(base[0] == w-1):
        dx=0
    dy = 1
    if(base[0] == h-1):
        dy = 0
    
    tl = [i for i in image.getpixel(tuple(base))]
    tlCoeff = (1.0-horCoeff) * (1.0-verCoeff)
    base[0] += dx
    tr = [i for i in image.getpixel(tuple(base))]
    trCoeff = horCoeff * (1.0-verCoeff)
    base[1] += dy
    br = [i for i in image.getpixel(tuple(base))]
    brCoeff = horCoeff * verCoeff
    base[0] -= dx
    bl = [i for i in image.getpixel(tuple(base))]
    blCoeff = (1.0-horCoeff) * verCoeff

    final = [  tr[i]*trCoeff + tl[i]*tlCoeff + br[i]*brCoeff + bl[i]*blCoeff for i in xrange(0, 3) ]
    return tuple(final)

# forces lat in range -pi/2 to pi/2, and lng to -pi to pi
# this is so buggy/inefficient but I don't want to figure out the 'right thing' right now.
def fix(lat, lng):
    while(lng < -math.pi):
        lng += 2.0 * math.pi
    while(lng > math.pi):
        lng -= 2.0 * math.pi
    return ( math.asin(math.sin(lat)), lng )


def placeSinuPixel(latp, lngp, w, h, image, p):
    (latp, lngp) = fix(latp, lngp)
    (xp, yp) = latLongToSinusoidal2d(latp, lngp, w, h)
    xp = int(round(xp))
    yp = int(round(yp))
    image.putpixel( (xp, yp), p )

# placeDisc()

# This is tricky to explain without pictures.
# There are two major concepts to get.

# 1) How I take a ball around lat/lng, and turn it into the guides.
# Basically, when two snake butts are sitting next to each other, I realized
# a ball of radius (path_width + field_width*.5) could be broken in half, each piece 
# translated to form one of the guides for a snake butt.  This is the basis of the
# cyanDots - they represent the centers of the original ball.

# 2) Elipses on a mercator project map to elipses on the surface of the sphere.  To get
# perfect circles on the sphere surface, I need to account for longitudinal skewing.
# That is the role of lngShift below.  One of the neatest things I learned at Harvard
# was that an ellipse can be thought of as treating the x coordinate as belonging to a circle
# with one radius, and the y coordinate as belonging to another circle with a different radius.
# As such, I basically just have to multiply r by lngShift when dealing with lng coordinates 
# to get the "right ellipse".

def placeDisc(x, y, image):
    (w,h) = image.size
    (lat, lng) = sinusoidal2dToLatLong(x, y, w, h)

    # I found these experimentally.  Really they should parameters
    # both are expressed in terms of degrees latitude along prime meridian
    pathWidth = .22
    fieldWidth = .08
    fieldR = fieldWidth * .5
    r = pathWidth + fieldR
    STEPS = 500
    for i in xrange(0, STEPS):
        theta = math.pi * 2.0 * i / STEPS 
        lngShift = 1.0 / math.cos(lat)
        latBase = lat
        lngBase = lng +  (pathWidth+fieldWidth)*lngShift*(1.0 if i > STEPS/2 else -1.0)
        latp = latBase + r * cos(theta)
        lngp = lngBase + r * sin(theta) * lngShift
        placeSinuPixel(latp, lngp, w, h, image, (255,255,255))
        # for a little extra help, mark the walls around the pivot too.
        latp = latBase + fieldR * cos(theta)
        lngp = lngBase + fieldR * sin(theta) * lngShift
        placeSinuPixel(latp, lngp, w, h, image, (255,255,255))

def transform(input_image, output_image):
    (w, h) = input_image.size
    assert (w,h) == output_image.size
    assert w == h*2
    dots = []
    for x in range(0, w):
        for y in range(0, h):
            p = input_image.getpixel( (x,y) )
            if ( p [0] < 10 and p[1] > 240 and p[2] > 240): #CYAN dot found
                dots.append( (x,y) )
            else:
                output_image.putpixel( (x,y), p )

    # We have to wait to place discs until full copy of image is done, or
    # the rest of the copy process above overwrites the disc placement.
    for d in dots:
        placeDisc(d[0], d[1], output_image)

def main ():
    assert len(sys.argv) == 3,\
        "Usage: %s <input_image> <output_image>" % sys.argv[0]

    # Get Ready
    input_image_file = sys.argv[1]
    output_image_file = sys.argv[2]

    # Read input image
    input_image = Image.open(input_image_file)
    if input_image.mode != 'RGB':
        input_image = input_image.convert("RGB")
    (width, height) = input_image.size
    assert width == height*2

    # Setup output image
    output_image = Image.new("RGB", (width, height), (255, 255, 255))

    # transform in -> out
    transform(input_image, output_image)

    # Finish
    output_image.save(output_image_file)

if __name__ == "__main__":
    main()
