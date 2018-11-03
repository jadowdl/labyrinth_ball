#########################################
# (c) 2011 Evan Mallory, TSM Games
#

import Image
import sys
from utils import *
from math import pi, cos

def main ():
	assert len(sys.argv) == 3,\
		"Usage: %s <output_image> <width>" % sys.argv[0]

	width = int(sys.argv[2])
	height = width / 2

	image = Image.new("RGB", (width, height), (70, 70, 70))

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
	image.save(sys.argv[1])

if __name__ == "__main__":
	main()
