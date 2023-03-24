"""
Image plane to generate the image from the traced 
photons trajectories in a curved spacetime

@author: AshCat
"""
from numpy import sqrt, sin, cos, arccos, arctan, pi, linspace, zeros
from common import Photon
from config import *

class image_plane:
    '''
    Given the initial cartesian coordinates in the image plane 
    (alpha,beta),the distance D to the force center and the 
    inclination angle i, this function calculates the initial 
    spherical coordinates (r, theta, phi) and the initial
    components of the momentum (kt, kr, ktheta, kphi)
    '''
    def __init__(self, D=D, iota = iota, s_side = screen_side, n_pixels = n_pixels):
        '''
        Defines a square NxN pixels screen with side of size s_side,
        located at a distance D and with an inclination iota.
        '''
        
        self.D = D 
        self.iota = iota  
        if n_pixels & 1:
            self.numPixels = n_pixels + 1
        else:
            self.numPixels = n_pixels 

        self.alphaRange = linspace(-s_side, s_side, self.numPixels)
        self.betaRange = linspace(-s_side, s_side, self.numPixels)
        print ("Size of the screen in Pixels:", self.numPixels, "X", self.numPixels)
        print ("Number of Pixels: ", (self.numPixels)**2)

    def photon_coords(self, alpha, beta, freq=1): 
        # Transformation from (Alpha, Beta, D) to (r, theta, phi) 
        r = sqrt(alpha**2 + beta**2 + self.D**2)
        theta = arccos((beta*sin(self.iota) + self.D*cos(self.iota))/r)
        phi = arctan(alpha/(self.D*sin(self.iota) - beta*cos(self.iota)))

        #Initial position of the photon in spherical coordinates 
        #(t=0, r, theta, phi)
        xin = [0., r, theta, phi]
                    
        
        #Given a frequency value w0, this calculates the initial 
        # 4-momentum of the photon  
        w0 =  freq    
        aux = alpha**2 + (-beta*cos(self.iota) + self.D*sin(self.iota))**2 
        kr =  (self.D/r)*w0   
        ktheta = (w0/sqrt(aux))*(-cos(self.iota) 
                    + (beta*sin(self.iota) + self.D*cos(self.iota))*(self.D/(r**2)))  
        kphi = - alpha*sin(self.iota)*w0/aux     
        kt = sqrt(kr**2 + r**2 * ktheta**2 + r**2*(sin(theta))**2 *kphi**2)
            
        # Initial 4-momentum in spherical coordinates  
        # (kt, kr, ktheta, kphi)
        kin = [kt, kr, ktheta, kphi]
        return xin, kin
    
    def create_photons(self):
        self.photon_list = []
        i=0
        for a in self.alphaRange:
            j = 0
            for b in self.betaRange:
                p = Photon(alpha=a, beta=b)
                p.xin, p.kin = self.photon_coords(a, b) 
                p.i, p.j = i, j
                p.initial_conditions()
                self.photon_list.append(p)
                j += 1
            i += 1
        return self.photon_list
    
    def create_image(self):
        self.image_data = zeros([self.numPixels,self.numPixels])
        for p in self.photon_list:
            self.image_data[p.i, p.j] = p.fP[1]

