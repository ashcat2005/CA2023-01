"""
===============================================================================
Schwarzschild metric

ds^2 = -(1-2M/r)dt^2 + dr^2 /(1-2M/r) + r^2 dtheta^2 + r^2 sin^2 (theta) dphi^2

===============================================================================
- Event horizon at r=2M
- ISCO at r = 6M

@author: Eduard Larrañaga - 2023
===============================================================================
"""

from numpy import sin, cos

class BlackHole:
    '''
    Definition of the Black Hole described by Schwarzschild metric
    '''
    def __init__(self, M):
        self.M = M
        self.EH = 2*M

    def metric(self,x):
        '''
        This procedure contains the Schwarzschild metric non-zero components in 
        spherical coordinates
        ===========================================================================
        Coordinates 
        t = x[0]
        r = x[1]
        theta = x[2]
        phi = x[3]
        ===========================================================================
        '''
        # Metric components
        gtt = -(1. - 2.*self.M/x[1])
        grr = 1./(1.- 2.*self.M/x[1])
        gthth = x[1]**2
        gphph = (x[1]*sin(x[2]))**2
        
        return [gtt, grr, gthth, gphph]

    def geodesics(self, q, lmbda):
        '''
        This function contains the geodesic equations in Hamiltonian form for 
        the Schwarzschild metric
        ===========================================================================
        Coordinates and momentum components
        t = q[0]
        r = q[1]
        theta = q[2]
        phi = q[3]
        k_t = q[4]
        k_r = q[5]
        k_th = q[6]
        k_phi = q[7]
        ===========================================================================
        Conserved Quantities
        E = - k_t
        L = k_phi
        ===========================================================================
        '''
        # Geodesics differential equations 
        dtdlmbda = q[4]*q[1]**2/(q[1]**2 - 2*self.M*q[1])
        drdlmbda = (1 - 2*self.M/q[1])*q[5]
        dthdlmbda = q[6]/q[1]**2
        dphidlmbda = q[7]/((q[1]*sin(q[2]))**2)
        
        dk_tdlmbda = 0.
        dk_rdlmbda = -self.M*(q[5]/q[1])**2 + q[6]**2/q[1]**3  \
                +q[7]**2/((q[1]**3)*sin(q[2])**2) \
                -self.M*(q[4]/(q[1]-2.*self.M))**2 
        dk_thdlmbda = (cos(q[2])/sin(q[2])**3)*(q[7]/q[1])**2
        dk_phidlmbda = 0.
        
        return [dtdlmbda, drdlmbda, dthdlmbda, dphidlmbda, 
                dk_tdlmbda, dk_rdlmbda, dk_thdlmbda, dk_phidlmbda]





###############################################################################

if __name__ == '__main__':
    print('')
    print('THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.')
    print('YOU NEED TO RUN THE main.py FILE TO GENERATE THE IMAGE')
    print('')