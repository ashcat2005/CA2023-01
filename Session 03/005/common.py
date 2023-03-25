"""
Common functions for ray tracing in a curved spacetime

@author: AshCat
"""
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