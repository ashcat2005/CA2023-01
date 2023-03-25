import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

from schwarzschild import *
from image_plane import *
from config import *


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
        color = 'crimson'
        sol = odeint(geodesics, p.iC, lmbda, args=(M,))
        indx = len(sol[:,1])
        for i in range(indx):
            if sol[i,1]<2.*M +1e-5: 
                indx = i
                color='black'
                break
        # Cartesian coordinates
        x = sol[:indx,1]*np.sin(sol[:indx,2])*np.cos(sol[:indx,3])
        y = sol[:indx,1]*np.sin(sol[:indx,2])*np.sin(sol[:indx,3])
        z = sol[:indx,1]*np.cos(sol[:indx,2])
        ax.plot(x, y, z, color=color)

    ax.set_xlim(-20,20)
    ax.set_ylim(-20,20)
    ax.set_zlim(-15,15)
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
    ax.set_zlabel(r'$z$')
    plt.show()




###### MAIN ######

detector = image_plane()

# Photons creation
photon_list = detector.create_photons()

# Integration and plot
plot_many3D(photon_list)
