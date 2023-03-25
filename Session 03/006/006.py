import matplotlib.pyplot as plt

from image_plane import *


def plot(image_data):
    '''
    Plots the image of the BH 
    '''
    ax = plt.figure().add_subplot(aspect='equal')
    ax.imshow(image_data.T, cmap = 'inferno', origin='lower')
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\beta$')
    plt.savefig('BH.jpeg')
    plt.show()


###### MAIN ######

detector = image_plane()

# Photons creation
detector.create_photons()

# Integration and plot
detector.create_image()
plot(detector.image_data)
