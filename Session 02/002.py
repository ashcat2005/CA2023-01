import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

from schwarzschild import *
from common import *
from image_plane import *


class Photon:
    def __init__(self, alpha, beta, detector, freq=1.):
        '''
        Given the initial coordinates in the image plane (alpha,beta), 
        this calculates the initial coordinates in spherical 
        coordinates (r, theta, phi)
        ''' 
        # Initial Cartesian Coordinates in the Image Plane
        self.alpha = alpha
        self.beta = beta

        #Initial position and momentum in spherical coordinates 
        self.xin, self.kin = detector.photon_coords(self.alpha, self.beta, freq)

        #Stores the initial conditions of coordinates and momentum 
        # to solve the geodesic equations.
        self.iC = None
        


def plot_many3D(photon_list):
    '''
    Plots the trajectory of the photon 
    '''
    # Independent parameter range
    lmbda = np.linspace(0,-150,1000)   

    ax = plt.figure().add_subplot(projection='3d')
    # Draw the black hole
    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    xs = 2*M*np.cos(u)*np.sin(v)
    ys = 2*M*np.sin(u)*np.sin(v)
    zs = 2*M*np.cos(v)
    ax.plot_surface(xs, ys, zs, color='k')

    for p in photon_list:
        sol = odeint(geodesics, p.iC, lmbda, args=(M,))
        indx = len(sol[:,1])
        for i in range(indx):
            if sol[i,1]<2.*M +1e-5: 
                indx = i
                break
        # Cartesian coordinates
        x = sol[:indx,1]*np.sin(sol[:indx,2])*np.cos(sol[:indx,3])
        y = sol[:indx,1]*np.sin(sol[:indx,2])*np.sin(sol[:indx,3])
        z = sol[:indx,1]*np.cos(sol[:indx,2])
        ax.plot(x, y, z, color='crimson')

    ax.set_xlim(-20,20)
    ax.set_ylim(-20,20)
    ax.set_zlim(-15,15)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$z$')
    plt.show()




###### MAIN ######

# Black Hole parameters
M = 1

# Detector  and creation of the detector object
D = 100*M
iota = np.pi/4
detector = image_plane(D=D, iota=iota)

# Photons definition
photon_list = []
# Initial Conditions in the image plane
ab_coords = [[ 5.,  5.],
             [-5., -5.],
             [ 5., -5.],
             [-5.,  5.],
             [0,0],
             [3,4]]
# Create photons
for i in range(len(ab_coords)):
    alpha, beta = ab_coords[i]
    p = Photon(alpha=alpha, beta=beta, detector=detector)
    p.iC = initCond(p.xin, p.kin, g, M)
    photon_list.append(p)

# Integration and plot
plot_many3D(photon_list)
