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
    lmbda = np.linspace(0,-150,1000)
    for p in photon_list:
        sol = odeint(geodesics, p.iC, lmbda, args=(M,))
        indx = len(sol[:,1])
        p.fP = [0,0,0,0,0,0,0,0]
        for i in range(indx):
            if sol[i,1]<2.*M +1e-5: 
                indx = i
                break
            elif sol[i,1]*cos(sol[i,2]) < 1e-5:
                indx = i
                p.fP = sol[i]
                break
    


###### MAIN ######

detector = image_plane()

# Photons creation
detector.create_photons()

# Integration and plot
geo_integ(detector.photon_list)
detector.create_image()
plot(detector.image_data)
