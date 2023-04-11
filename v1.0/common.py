"""
===============================================================================
Common functions for ray tracing in a curved spacetime
===============================================================================
@author: Eduard Larrañga - 2023
===============================================================================
"""
from scipy.integrate import odeint
from numpy import linspace, cos, zeros
import matplotlib.pyplot as plt
import sys


def initCond(x, k, blackhole):
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
    g_tt, g_rr, g_thth, g_phph = blackhole.metric(x)
    
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
    
    def initial_conditions(self, blackhole):
        '''
        Calcualtes and stores the initial conditions of coordinates 
        and momentum to solve the geodesic equations.
        '''
        self.iC = initCond(self.xin, self.kin, blackhole)


def geo_integ(p, blackhole, in_edge, out_edge):
    '''
    Integrates the motion equations of the photon 
    '''
    lmbda = linspace(0,-200,1500)
    sol = odeint(blackhole.geodesics, p.iC, lmbda)
    indx = len(sol[:,1])
    p.fP = [0,0,0,0,0,0,0,0]
    for i in range(indx):
        if sol[i,1] < blackhole.EH + 1e-5: 
            indx = i
            break
        elif abs(sol[i,1]*cos(sol[i,2])) < 1e-1:
            if sol[i,1] > in_edge and sol[i,1] < out_edge:
                indx = i
                p.fP = sol[i]
                break

class Image:
    '''
    Image class
    Creates the photon list and generates the image
    '''
    def __init__(self):
        pass

    def create_photons(self, blackhole, detector):
        '''
        Creates the photon list
        '''
        self.detector = detector
        self.photon_list = []
        i=0
        for a in self.detector.alphaRange:
            j = 0
            for b in self.detector.betaRange:
                p = Photon(alpha=a, beta=b)
                p.xin, p.kin = self.detector.photon_coords(a, b) 
                p.i, p.j = i, j
                p.initial_conditions(blackhole)
                self.photon_list.append(p)
                j += 1
            i += 1
    
    def create_image(self, blackhole, acc_structure):
        '''
        Creates the image data 
        '''
        self.image_data = zeros([self.detector.numPixels, self.detector.numPixels])
        photon=1
        for p in self.photon_list:
            geo_integ(p, blackhole, acc_structure.in_edge, acc_structure.out_edge)
            self.image_data[p.i, p.j] = acc_structure.spectrum(p.fP[1])
            sys.stdout.write("\rPhoton # %d" %photon)
            sys.stdout.flush()
            #print('Photon #',photon)
            photon +=1

    def plot(self, savefig=False, filename=None):
        '''
        Plots the image of the BH 
        '''
        ax = plt.figure().add_subplot(aspect='equal')
        ax.imshow(self.image_data.T, cmap = 'inferno', origin='lower')
        ax.set_xlabel(r'$\alpha$')
        ax.set_ylabel(r'$\beta$')
        if savefig:
            plt.savefig(filename)
        plt.show()




###############################################################################

if __name__ == "__main__":
    print('')
    print('THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.')
    print('YOU NEED TO RUN THE main.py FILE TO GENERATE THE IMAGE')
    print('')