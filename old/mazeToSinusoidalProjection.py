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

# see discussion of "lightness" - http://en.wikipedia.org/wiki/HSL_and_HSV
def getBrightness(v):
    if type(v) is tuple:
        return (v[0]*.3 + v[1]*.59 + v[2]*.11)/255.0
    else:
        return float(v)/255.0

def mapsToField(lat, lng, input_image):
    (w,h) = input_image.size
    (x,y) = latLongToRegular2d(lat, lng, w, h)
    assert 0 <= x and x < w and 0 <= y and y < h, "X=%s, Y=%s"%(x,y)
    # print lat, lng, x, y
    valence = getBrightness(input_image.getpixel( (int(x), int(y)) ))
    return valence < .3

def overlayField(input_image, output_image):
    (w, h) = output_image.size
    for x in range(0, w):
        for y in range(0, h):
            (lat, lng) = sinusoidal2dToLatLong(x, y, w, h)
            if( fabs(lng) > pi): continue # given lat, -pi to pi corresponds
                                       # to those x in range on stencil.
            if(mapsToField(lat, lng, input_image)):
                output_image.putpixel((x,y), (140,140, 140)) # light gray

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
        input_image = image.convert("RGB")
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
