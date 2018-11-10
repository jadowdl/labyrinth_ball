#!/usr/bin/env bash

# Make sure you use venv...
if [ -d ".venv" ]
then
  source .venv/bin/activate
else
  echo "No .venv/ folder found, please initialize virtualenv for python2"
  echo "(e.g., 'pip install virtualenv && virtualenv -p /usr/bin/python2.7 .venv')"
  exit 1
fi

# Update pip dependencies
pip install -r pip_requirements.txt

# Uncomment the example you want to run.

# Hand-crafted 5 circuit labyrinth
python labyrinth_ball/interactive/makeModelFromSineMap.py labyrinth_ball/res/labyrinth_5.png

# Hand-crafted 7 circuit labyrinth
# python labyrinth_ball/interactive/makeModelFromSineMap.py labyrinth_ball/res/labyrinth_7.png

# Exit virtualenv
deactivate
