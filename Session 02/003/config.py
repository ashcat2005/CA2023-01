"""
Parameters to create the image

@author: AshCat
"""

from numpy import pi
from schwarzschild import *



# Black Hole parameters
M = 1

# Detector parameters
D = 100*M
iota = pi/4
screen_side = 5
n_pixels = 20

# Metric
metric = g
