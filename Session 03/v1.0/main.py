from image_plane import *
from common import plot


###### MAIN ######

detector = image_plane()

# Photons creation
detector.create_photons()

# Create the image data
detector.create_image()

# Plot the image
plot(detector.image_data, savefig=True)
