"""
Thind disk model

@author: AshCat
"""


def thin_disk(r, R_min , R_max):
    m = (1.-0.)/(R_min - R_max)
    intens = m * (r - R_max)
    if r>R_min and r<R_max:
        return intens
    else:
        return 0.
