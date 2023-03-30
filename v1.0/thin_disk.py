"""
===============================================================================
Thind accretion disk with a simple linear model of spectrum
===============================================================================
@author: Eduard Larrañaga - 2023
===============================================================================
"""

class thin_disk:
    def __init__(self, R_min , R_max):
        self.in_edge = R_min
        self.out_edge = R_max

    def spectrum(self, r):
        m = (1.-0.)/(self.in_edge - self.out_edge)
        intensity = m * (r - self.out_edge)
        if r>self.in_edge and r<self.out_edge:
            return intensity
        else:
            return 0.



###############################################################################

if __name__ == '__main__':
    print('')
    print('THIS IS A MODULE DEFINING ONLY A PART OF THE COMPLETE CODE.')
    print('YOU NEED TO RUN THE main.py FILE TO GENERATE THE IMAGE')
    print('')
