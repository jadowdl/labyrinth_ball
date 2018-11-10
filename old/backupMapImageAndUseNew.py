#########################################
# (c) 2011 Evan Mallory, TSM Games
#

# This script takes the name of a model sinusoidal map image, like
# "labyrinth_7.png" as an argument. It moves that file into a backup directory
# after adding a random postfix, then copies "out-<name>" (that is created
# by makeModelFromSineMap.py) to <name>.
#
# Note - Nov 9 2018, James Dowdell
#
# I've peeked at this and I don't think we need it, because fbofw we're storing
# images in git now and that will prevent clobbering.  Easy enough to recreate
# anyway.

import sys
import os
from utils import *
from math import sqrt
import random
import shutil

def randomString(size):
	result = ""
	for i in range(0, size):
		result += chr(random.randint(0, 25) + 65)
	return result

def main():
	random.seed()

	assert len(sys.argv) == 2,\
		"Usage: %s <map_image>" % sys.argv[0]

	old = sys.argv[1]
	new = "out-" + old

	if not os.path.exists("backup"):
		os.mkdir("backup")
	backup = (os.path.join("backup", os.path.splitext(old)[0]) + "_" +
		randomString(4) + ".png")

	shutil.move(old, backup)
	shutil.move(new, old)

if __name__ == "__main__":
	main()

