import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

from schwarzschild import *
from image_plane import *


def plot(image_data):
    '''
    Plots the image of the BH 
    '''
    
    ax = plt.figure().add_subplot(aspect='equal')
    ax.imshow(image_data, cmap = 'inferno')
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\beta$')
    plt.show()


def geo_integ(photon_list):
    '''
    Integrates the motion equations of the photon 
    '''
    # Independent parameter step
    d_lmbda = 0.05 

    photon =1
    for p in photon_list:
        print('photon ',photon)
        fP = p.iC
        i=0
        while fP[1]*cos(fP[2]) > 1e-5 and fP[1]>2.*M +1e-10 and i<100000:
            step = odeint(geodesics, fP, [0, -d_lmbda], args=(M,))
            fP = step[-1]
            i+=1
        if i>99999:
            fP = np.NaN
        p.fP = fP
        photon += 1


###### MAIN ######

detector = image_plane()

# Photons creation
detector.create_photons()

# Integration and plot
geo_integ(detector.photon_list)
detector.create_image()
plot(detector.image_data)
