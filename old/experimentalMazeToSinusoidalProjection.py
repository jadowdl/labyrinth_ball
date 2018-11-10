#########################################
# (c) 2011 James Dowdell, TSM Games
#

import Image
import sys
from utils import *
from math import pi, cos, fabs

def fillStencil(image):
    (width, height) = image.size 
    lat_count = 0
    long_count = 0
    for latitude in drange(-pi/2, pi/2 + .01, pi/32):
        lat_count += 1
        long_count = 0
        for longitude in drange(-pi, pi + .01, pi/32):
            long_count += 1
            (x, y) = latLongToSinusoidal2d(latitude, longitude, width, height)
            if lat_count % 4 == 1 or long_count % 4 == 1:
                image.putpixel((x, y), (256, 0, 0))
            else:
                grey_color = 128
                image.putpixel((x, y), (grey_color, grey_color, grey_color))

#def lookupAliasedPixel(x, y, image):
#    (w,h) = image.size
#    x *= float(w-1)/w
#    y *= float(h-1)/h
#    x = int(x)
#    y = int(y)
#    return image.getpixel( (x,y) )

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

def classifyPixelXY(x, y, image):
    v = lookupAliasedPixel(x, y, image)
    isProtected = v[1] > .7 # green channel is the protected channel
    bright = v[0] # getBrightness(v) TODO
    if(isProtected):
        bright = 255.0
    return [bright, isProtected]

def classifyPixel(r, theta, image):
    (w,h) = image.size
    assert r <= w*.5, "R and W %s %s"%(r, w)
    return classifyPixelXY( .5*w + r*cos(theta), .5*h + r*sin(theta), image)

def overlayOnSinusImage(lat, long, brightness, image, step=0.0):
    if(step > .005):
        k = int(step / .005)
        for i in xrange(0, k):
            overlayOnSinusImage(lat, long, brightness, image)
            long += .005
        return
    (x,y) = latLongToSinusoidal2d(lat, long, image.size[0], image.size[1])
    #print lat, long, x,y
    brightness = int(brightness)
    image.putpixel( (x,y), (brightness, brightness, brightness) )

def overlayField(input_image, output_image):
    (w, h) = input_image.size
    assert w == h
    for r in range(1, w): # for each circle, we'll project it to the output image
        radius = r*.5 # cheap trick; really want to step by .5 from 0 to width/2
        sigma_step = .5/radius # traveling sigma_step radians should be 
                               # less than .5 pixels.
        assert sigma_step > .000001 # want to finish the loop below
        # Part A - step through all sectors on the circle and classify them
        sectors = []
        theta = 0.0
        count = 0
        protectedCount = 0
        while(theta < 2.0 * pi ):
            t = classifyPixel(radius, theta, input_image) # (fieldValence, isProtected)
            # DEBUGGING t[1] = False  # TODO
            count += 1
            theta += sigma_step
            if(t[1]):
                protectedCount+=1
            sectors.append(t)
        # Part B - figure out how you're going to shrink protected and not protected sectors
        (lat, long) = regularRThetaToLatLong(radius, 0.0,w,h) # just want lat
        latMult = pi # d units in 2D space from center to extreme becomes pi units on the sphere from north to south pole.
        latCircumference = cos(lat) * 2.0 * pi
        prtMult = latMult / latCircumference # at equator, no shrink; towards poles, we've shrunk by latCircumference, so counter.

        # protectedSectorLength * numProtected + regularSectorLength * numNotProtected = latCircumference on unit sphere
        # => latMult * sigma_step * protectedCount + mu * sigma_step*(count-protectedCount) = latCircumference
        mu = 0.0
        print count, protectedCount
        if(protectedCount != count):
            #mu = ((latCircumference / sigma_step) - prtMult * protectedCount) / (count-protectedCount) # TODO
            #mu = (2.0*pi/sigma_step - prtMult*protectedCount)/(count-protectedCount)
            mu = (latCircumference / (sigma_step * pi) - protectedCount) / (count-protectedCount)

        # Part C - shrink/stretch the sectors and map to the output_image
        # protected sectors stretch by latMult, so that squares map to squares
        # non protected sectors stretch by mu
        long = -pi 
        for t in sectors:
            step = 0.0
            if(t[1]):
                step = 2.0*pi/latCircumference    * (sigma_step * pi)
            else:
                step = 2.0*pi/latCircumference    * (sigma_step * pi * mu)
            overlayOnSinusImage(lat, long, t[0], output_image, step)
            long += step

def main ():
    assert len(sys.argv) == 4,\
        "Usage: %s <input_image> <output_image> <width>" % sys.argv[0]

    # Get Ready
    input_image_file = sys.argv[1]
    output_image_file = sys.argv[2]
    width = int(sys.argv[3])
    assert width%2==0
    height = width / 2
    output_image = Image.new("RGB", (width, height), (70, 70, 70))

    # Read input image
    input_image = Image.open(input_image_file)
    if input_image.mode != 'RGB':
        input_image = input_image.convert("RGB")
    (width, height) = input_image.size
    assert width == height

    # fillStencil
    fillStencil(output_image)

    # Overlay field
    overlayField(input_image, output_image)

    # Finish
    output_image.save(output_image_file)

if __name__ == "__main__":
    main()
