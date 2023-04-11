"""
===============================================================================
Main script 
Creates the Black Hole image
===============================================================================
@author: Eduard Larrañga - 2023
===============================================================================
"""

from common import Image
from config import *



#################################### MAIN #####################################

image = Image()

# Photons creation
image.create_photons(blackhole, detector)

# Create the image data
image.create_image(blackhole, acc_structure)

# Plot the image
image.plot(savefig=True, filename=filename)
