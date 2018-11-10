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


# see discussion of "lightness" - http://en.wikipedia.org/wiki/HSL_and_HSV
def getBrightness(v):
    if type(v) is tuple:
        return (v[0]*.3 + v[1]*.59 + v[2]*.11)/255.0
    else:
        return float(v)/255.0

def transform(input_image, output_image):
    (w, h) = input_image.size
    assert (w,h) == output_image.size
    assert w == h
    for r in range(1, w): # for each circle, we'll project it to the output image
        radius = r*.5 # cheap trick; really want to step by .5 from 0 to width/2
        sigma_step = .5/radius # traveling sigma_step radians should be 
                               # less than .5 pixels.
        assert sigma_step > .000001 # want to finish the loop below
        
        # setup for transformation on this circle
        (lat, long) = regularRThetaToLatLong(radius, 0.0,w,h) # just want lat
        adj = radius * math.pi / ( cos(lat) *  w* .5)
        if(adj > 7.0):
            print r, "Too much distortion requested"
            adj = 7.0

        # step through all sectors on the circle
        theta = 0.0
        print r
        while(theta < 2.0 * math.pi ):
            # pixel in 
            x = w*.5 + radius * cos(theta)
            y = h*.5 + radius * sin(theta)
            pixel = lookupAliasedPixel(x,y,input_image)
            pixel = tuple( [int(pixel[i]) for i in xrange(0, 3)] )

            # transform (r, theta) to (r', theta')
            ######################################

            # treat as t between 0 and pi.
            t = theta
            neg = False
            if(t > math.pi):
                t = math.fabs(t-2.0*math.pi)
                neg = True
            # piece wise squishing.
            if(t < math.pi/14):
                t *= adj
            elif(t > 13.0 * math.pi/14): # with north skew, have to protect western sector!
                t = math.pi - (math.pi-t) * adj
            else:
                # t' = pi - (pi-c)/(pi-b) * (pi-t)
                t  = .5*math.pi - (.5*math.pi - adj * math.pi/14.0) / (.5*math.pi - math.pi/14.0) * (.5*math.pi - t)
            # undo negation if any.
            if(neg):
                t *= -1
            new_r = radius
            new_theta = t

            # pixel out
            x = w*.5 + new_r * cos(new_theta)
            y = h*.5 + new_r * sin(new_theta)
            x *= float(w-1)/w
            y *= float(h-1)/h
            output_image.putpixel( (int(round(x)),int(round(y))), pixel)

            #
            theta += sigma_step

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
    assert width == height

    # Setup output image
    output_image = Image.new("RGB", (width, height), (255, 255, 255))

    # transform in -> out
    transform(input_image, output_image)

    # Finish
    output_image.save(output_image_file)

if __name__ == "__main__":
    main()
