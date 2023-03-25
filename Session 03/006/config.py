"""
Parameters to create the image

@author: AshCat
"""

from numpy import pi
from schwarzschild import *
from thin_disk import *


# Metric
metric = g
geodesics = geodesics

# Black Hole parameters
M = 1

# Detector parameters
D = 100*M
iota = pi/2.5
screen_side = 20*M
n_pixels = 100

# Disk Model
disk = thin_disk
R_min = 6*M
R_max = 13*M
