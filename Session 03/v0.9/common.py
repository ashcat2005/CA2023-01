"""
Common functions for ray tracing in a curved spacetime

@author: AshCat
"""
from scipy.integrate import odeint
from numpy import linspace
import matplotlib.pyplot as plt

from config import *

def initCond(x, k, metric=metric, M=M):
    '''
    Given the initial conditions (x,k)
    this function returns the list
    [t, r, theta, phi, k_t, k_r, k_theta, k_phi] 
    with the initial conditions needed to solve 
    the geodesic equations 
    (with the covariant components of the momentum vector)    
    # Coordinates
    t = x[0]
    r = x[1]
    theta = x[2]
    phi = x[3]
    kt = k[0]
    kr = k[1]
    ktheta = k[2]
    kphi = k[3]
    '''
    # Metric components
    g_tt, g_rr, g_thth, g_phph = metric(x, M=M)
    
    # Lower k-indices
    k_t = g_tt*k[0] 
    k_r = g_rr*k[1]
    k_th = g_thth*k[2]
    k_phi = g_phph*k[3]
    
    return [x[0], x[1], x[2], x[3], k_t, k_r, k_th, k_phi]



class Photon:
    def __init__(self, alpha, beta, freq=1.):
        '''
        Given the initial coordinates in the image plane (X,Y), the distance D 
        to the force center and inclination angle i, this calculates the 
        initial coordinates in spherical coordinates (r, theta, phi)
        ''' 
        # Initial Cartesian Coordinates in the Image Plane
        self.alpha = alpha
        self.beta = beta
        
        # Pixel coordinates
        self.i = None
        self.j = None

        #Initial position and momentum in spherical coordinates 
        self.xin = None 
        self.kin = None 
        
        # Stores the final values of coordinates and momentum 
        self.fP = None
    
    def initial_conditions(self):
        '''
        Calcualtes and stores the initial conditions of coordinates 
        and momentum to solve the geodesic equations.
        '''
        self.iC = initCond(self.xin, self.kin)


def geo_integ(p):
    '''
    Integrates the motion equations of the photon 
    '''
    lmbda = linspace(0,-150,1000)
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

def plot(image_data, savefig=False):
    '''
    Plots the image of the BH 
    '''
    ax = plt.figure().add_subplot(aspect='equal')
    ax.imshow(image_data.T, cmap = 'inferno', origin='lower')
    ax.set_xlabel(r'$\alpha$')
    ax.set_ylabel(r'$\beta$')
    if savefig:
        plt.savefig(filename)
    plt.show()